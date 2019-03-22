function [fgMaskMov,bgvals] = get_foreground_mask_movie(movieRotCyc, foregroundMaskingSettings)
    % get_foreground_mask_movie
    
    % :param movieRotCyc: rotated movie array
    % :param foregroundMaskingSettings: settings
    %
    % :returns: fgMaskMov
    
    % this whole procedure should be made more robust by choosing a better
    % method for selecting signal/background threshold. Right now 
    % first threshold- graytresh on all the non-nan values
    % second threshold - mean(background)+3std (maybe more?)
    % how to deal with anomalies?
    % Should check how this works for regular data (i.e. how many errors
    % etc)
    
    import AB.Processing.amplify_movie;

        
    % simplify parameter name
    maxAmpDist = foregroundMaskingSettings.maxAmpDist;
    
    % compute a nan mask, rather could pass as parameter since computed
    % before
    nanMask = isnan(movieRotCyc);
    
    % signal pixels
    movMask = ~nanMask;
    
    % statistics, i.e. mean and std. Todo: potential failure if there is a
    % lot of molecules and not a lot of background, since the mean is taken
    % over everything
    movMeanVal = mean(movieRotCyc(movMask));
    movMeanStd = std(movieRotCyc(movMask));
    
    % num values that are non-defined
    numRandVals = sum(nanMask(:));
    
    currRng = rng(); % so that current rng state can be restored
    % temporarily set to produce predictable pseudorandom values for reproducibility
    rng(rng(0, 'twister'));  
    randVals = randn([numRandVals, 1]) .* movMeanStd + movMeanVal;
    rng(currRng); % restore rng state
    
    % assign random values to the nan's, so we could work with the whole
    % image directly
    movieRotCycAmp = movieRotCyc;
    movieRotCycAmp(nanMask) = randVals;
    
    szInit = size(movieRotCycAmp);

    % amplify movie using maxAmpDist.
    % movie size changes because of amplification
    [movieRotCycAmp] = amplify_movie(movieRotCycAmp, maxAmpDist);
    
    % there is going to be a invalid region because the filters don't use
    % the values beyond the edges. We should already know this from the
    % filter size, this assumes that they are padded on both sides equally
    invalidPadSz = (szInit - size(movieRotCycAmp))./2; % since movieRotCycAmp is restricted to valid region
    
    % substract the mean from the amplified movie
    a = movieRotCycAmp - mean(movieRotCycAmp(:));
    
    % Take the sum of the values, leave only x direction
	colSignalVect = sum(sum(a,3), 1); 
    
    % these should come from input parameters
    maxSigmaNonBlip = foregroundMaskingSettings.maxSigmaNonBlip;
    %minValDistBetweenAdjLocalExtrema = 0;
    
    % ignore the rows having negative values
    colSignalVect(colSignalVect<0) = nan;
    
    % find peaks, use inbuild findpeaks of matlab
    [PKS,LOCS,W,P] = findpeaks(colSignalVect);
    
    % define signal channels mask
    sigChannelsMask = zeros(1,length(colSignalVect));
    
    % if first location is closer than the edge, remove it
    if LOCS(1) < maxSigmaNonBlip 
        LOCS = LOCS(2:end);
    end
    
    % the same for last location
	if LOCS(end) > (length(colSignalVect)-maxSigmaNonBlip+1)
        LOCS = LOCS(1:end-1);
    end
    
    % assuming width of the channel is 2*maxSigmaNonBlip+1, we bitmask the
    % channel mask
    for i=1:length(LOCS)
        sigChannelsMask(max(1,LOCS(i)-maxSigmaNonBlip):min(LOCS(i)+maxSigmaNonBlip,length(sigChannelsMask)))= 1;
    end
    
    % pad x dimension with 0's at left and right
	sigChannelsMask = padarray(sigChannelsMask, [0, invalidPadSz(2)], false, 'both');

    % pad amplified movie with nan's
	nonnanSigChannelVals = padarray(movieRotCycAmp, invalidPadSz, NaN, 'both');
    % also the places previously assigned values, assign nan's again
    nonnanSigChannelVals(nanMask) = NaN;
    
    allValues = nonnanSigChannelVals(not(isnan(nonnanSigChannelVals )));
   fgThreshA1 = graythresh(allValues); %TODO: make this potential point of failure more robust

    % select values in the channels ?? why, but in the channel we expect to
    % be only signal
    nonnanSigChannelVals = nonnanSigChannelVals(:, not(logical(sigChannelsMask)), :);
    nonnanSigChannelVals = nonnanSigChannelVals(not(isnan(nonnanSigChannelVals)));
    fgThreshA2 = mean(nonnanSigChannelVals)+3*std(nonnanSigChannelVals);
    
    bgvals = [mean(nonnanSigChannelVals) std(nonnanSigChannelVals)];
    % threshold inside the channels
%     fgThreshA1 = graythresh(nonnanSigChannelVals); %TODO: make this potential point of failure more robust
    
    %
%     bgValsA1 = nonnanSigChannelVals(nonnanSigChannelVals < fgThreshA1);
%     
%     % second threshold
%     fgThreshA2 = mean(bgValsA1) + 3*std(bgValsA1);
%     
    % should we take max or min?
    fgThreshB = max(fgThreshA1, fgThreshA2);
    fgMaskMov =  padarray(movieRotCycAmp >= fgThreshB, invalidPadSz, false, 'both');
    
   
end