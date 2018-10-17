% function featureDistances=getFeatureDistsances(inputStruct_smooth,inputStruct_raw)
function [featuresCellArray_processed, featureDistances, featureDistanceVariances] = get_feature_distances(kymo,settings,backgroundCompensation)
%     import dotkymoAlignment.*
%     settings=settings_dotKymoAlignment; % A struct containing all parameters used.
    import ELD.Core.get_color_map;

    if nargin < 1 || isempty(kymo)
        
        import ELD.Import.import_tiff_img
        kymo = import_tiff_img();
        
%         [kymoFilenames, dirpath] = uigetfile('*.mat', 'Select mat-file with unaligned kymograph(s)', 'Multiselect', 'on');
%         aborted = isequal(dirpath, 0);
%         if aborted
%             return;
%         end
%     %     if not(iscell(kymoFilenames))
%     %         kymoFilenames = {kymoFilenames};
%     %     end
%         kymoFilepaths = fullfile(dirpath, kymoFilenames);
% 
%         kymo = load(kymoFilepaths);
%         kymo = kymo.processedKymo;
%         
%         kymo = arrayfun(@(input) input.kymo_dynamicMeanSubtraction,kymo,'uniformoutput',false);
%     %     if nargin < 1
%     %         inputStructs = load(kymoFilepaths);
%     %     end
% 
%     % elseif size(inputStruct_smooth) ~= size(inputStruct_raw)
%     %     return;
%     else if isstruct(kymo)
%         kymo = kymo.processedKymo;
%         kymo = arrayfun(@(input) input.kymo_dynamicMeanSubtraction,kymo,'uniformoutput',false);
%         else
%             kymo = kymo;
%         end
    end
%     
%     if nargin < 2 || isempty(minOverlap)
%         minOverlap = 5;
%     end
%     
%     if nargin < 3 || isempty(confidenceInterval)
%         confidenceInterval = 1.5;
%     end
    
    if nargin < 2 || isempty(settings)
        import ELD.Import.load_eld_kymo_align_settings;
        settings = load_eld_kymo_align_settings();
    end
    
    if nargin < 3 || isempty(backgroundCompensation)
        backgroundCompensation = false;
    end
    
    minOverlap = settings.minVerticalOverlap;
    confidenceInterval = settings.confidenceInterval;
    
    
    
%     if backgroundCompensation
%         import ELD.Core.generate_smart_kymo_masks;
%         [ ~ , backgroundBitmask ] = generate_smart_kymo_masks( kymo, 1, settings.maxGapLen, settings.minSegmentLen, settings.dilationFactor, settings.numThresholds );

%         signalImg = kymoImg.*signalBitmask;
% 
%         signalImg(signalImg == 0) = NaN;

        % [rowInd,colInd] = find(quantImg);
        % signalImg = signalImg(:,min(colInd):max(colInd));

%         figure(1)
%         subplot(2,2,4);
%         imagesc(signalImg);

        % s=signalImg(signalImg>0)
        % barcode = nanmean(signalImg(signalImg>0))
        % barcode = nansum(signalImg)./nansum(signalImg~=0);
%         barcode = nanmean(signalImg);
        % quantImg == 0
        % qm = kymoImg.*(quantImg == 0);
%         backgroundImg = kymo.*backgroundBitmask;
%         backgroundImg(backgroundImg == 0) = NaN;
%         backgroundFramewise = nanmean(backgroundImg,2);
%     end

%     for row = 1:size(kymo,1)
%         kymo(row,:) = kymo(row,:) - mean(kymo(row,:));
%     end
%     kymo = kymo{1};
%     kymo = kymo - mean(kymo(:));
%     kymo = kymo / max(kymo(:));
    
%     figure, imagesc(kymo);


    %After the background subtraction, a kymogrpah may have negative values which 
    %can't be accurately displayed so a corresponding kymograph without any
    %background subtraction
    %is used for all visual demonstration purposes. All such results are marked
    %with the suffix "_show"

    % raw=inputStruct_raw.kymo_dynamicMeanSubtraction;
    % raw_show=inputStruct_raw.kymo_noSubtraction;
    
%     numKymos = length(kymo);
%     numKymos = 1;
%     kymosToRun = [3];

    fullTimer = tic;
    
%     kymosToRun = 1:numKymos;
%     for kymoNum = kymosToRun
        
%         kymoTimer = tic;
        
%         peakTimer1 = tic;
%         imgArr=inputStruct_smooth(kymoNum).kymo_dynamicMeanSubtraction;
%         imgArr = kymoImgArrs{kymoNum};
%         imgArr_show=inputStruct_smooth(kymoNum).kymo_noSubtraction;

    import ELD.Core.generate_peak_bitMap;
    import ELD.Core.generate_peak_graph;
    import ELD.Core.process_peak_graph;
    import ELD.Core.organize_features;
    if backgroundCompensation
        peakBitMap=generate_peak_bitMap(kymo);
        peakBitMap=peakBitMap.*(~(kymo<=0));


        %Remove all isolated pixels in the peakBitMap i.e pixles which dont have pixel above or below within the range specified by 'k'.
        sel=zeros(3,2*settings.localFluctuationWindow+1);
        sel(1,:)=1;
        sel(3,:)=1;
        tmp_1=imdilate(peakBitMap,sel); 
        peakBitMap=peakBitMap.*tmp_1;
%         peakTime1(kymoNum) = toc(peakTimer1)

%         graphTimer1 = tic;
        %Graph Based Computatations.

        S=generate_peak_graph(kymo,settings.localFluctuationWindow,peakBitMap);%Create a graph
%         graphTime1(kymoNum) = toc(graphTimer1)
%         processingTimer1 = tic;
        %Process the connected components and save only those whose size is greater than the threshold.

        featuresCellArray=process_peak_graph(S,settings.minimumSizeOfConnectedComponent);
%         processingTime1(kymoNum) = toc(processingTimer1)
%         orderTimer1 = tic;

        if size(featuresCellArray,2) < 2
            fprintf('Less than two features were detected.\n');
            featuresCellArray_processed = [];
            featureDistances = [];
            featureDistanceVariances = [];
            return;
        end

        %Organize components by mean horizontal position, and generate a
        %feature map for further computations, and a color map for displaying purposes.

        [featureMap,colorMap,featuresCellArray_ordered] = organize_features(featuresCellArray,size(kymo));
%         orderTime1(kymoNum) = toc(orderTimer1)
        
        
        %%Backgroundsubtraction

%         backgroundTimer = tic; 
        signalMap = featureMap >= 1;
        signalMap = imdilate(signalMap,ones(1,11));

        import ELD.Core.subtract_kymo_background
        kymo = subtract_kymo_background( kymo, signalMap );  
    end
        
%         backgroundNoise = nanmean(imgArr.*~signalMap);
%         figure, plot(backgroundNoise);
% 
%         coeff = polyfit(1:length(backgroundNoise),backgroundNoise,1);
% 
%         backgroundNoise = coeff(2) + coeff(1) .* (1:length(backgroundNoise));
%         figure, plot(backgroundNoise);
% 
%         imgArr_backGroundSubtracted = imgArr;
%         for row = 1:size(imgArr,1)
%             imgArr_backGroundSubtracted(row,:) = imgArr(row,:) - backgroundNoise;
%         end
    
%         backgroundTime(kymoNum) = toc(backgroundTimer)

%         peakTimer2 = tic; 

    peakBitMap=generate_peak_bitMap(kymo);
    peakBitMap=peakBitMap.*(~(kymo<=0));


    %Remove all isolated pixels in the peakBitMap i.e pixles which dont have pixel above or below within the range specified by 'k'.
    sel=zeros(3,2*settings.localFluctuationWindow+1);
    sel(1,:)=1;
    sel(3,:)=1;
    tmp_1=imdilate(peakBitMap,sel); 
    peakBitMap=peakBitMap.*tmp_1;
%         peakTime2(kymoNum) = toc(peakTimer2)

%         graphTimer2 = tic;
    %Graph Based Computatations.
    S=generate_peak_graph(kymo,settings.localFluctuationWindow,peakBitMap);%Create a graph
%         graphTime2(kymoNum) = toc(graphTimer2)
%         processingTimer2 = tic;
    featuresCellArray=process_peak_graph(S,settings.minimumSizeOfConnectedComponent);%Process the connected components and save only those whose size is greater than the threshold.

    if size(featuresCellArray,2) < 2
        fprintf('Less than two features were detected.\n');
        featuresCellArray_processed = [];
        featureDistances = [];
        featureDistanceVariances = [];
        return;
    end
    %         processingTime2(kymoNum) = toc(processingTimer2)
%         orderTimer2 = tic;
    [featureMap,colorMap,featuresCellArray_ordered] = organize_features(featuresCellArray,size(kymo));
%         orderTime2(kymoNum) = toc(orderTimer2)






%         figure;
%         sp(1) = subplot(1,2,1);
% %         imagesc(imgArr_backgroundSubtracted);
%         imagesc(imgArr);
%         colorbar;
%         sp(2) = subplot(1,2,2);
%         imagesc(colorMap);
%         
%         linkaxes(sp,'xy');


    numFeatures = numel(featuresCellArray_ordered);
%         numColors = numFeatures;
%         twoFeatureOverlap = zeros(numFeatures);
%         featureDistance = zeros(numFeatures);
%         featureSquareDistance = zeros(numFeatures);

    distancesToCheck = 1:numFeatures;

%         tempTimer = tic;
%         [featureDistances, featureDistanceVariances , featureOverlaps] = calculate_feature_distances(featuresCellArray_ordered , 9 );

    import ELD.Core.calculate_feature_distances
    [featureDistance, featureDistanceVariance, featureOverlaps] = calculate_feature_distances(featuresCellArray_ordered,distancesToCheck,[],minOverlap);
%         tempTime = toc(tempTimer)

    import ELD.Core.merge_features;
    import ELD.Core.order_features_by_dists;
%     round = 1;

    while true


%         disp(num2str(round));
%         round = round+1;
%             estimatorTimer = tic;

        optimalEstimatorMembers = cell(numFeatures);
        optimalEstimatorDistance = nan(numFeatures);
        optimalEstimatorVariance = nan(numFeatures);
        displayCells = cell(numFeatures);
        for featureI = 1:numFeatures-1
            for featureJ = featureI+1:min(numFeatures,featureI+5)
                [optimalEstimatorVariance(featureI,featureJ), optimalEstimatorMembers{featureI,featureJ}, ~] = graphshortestpath(sparse(featureDistanceVariance),featureI,featureJ);

                for featureK = 1:length(optimalEstimatorMembers{featureI,featureJ})-1
%                         optimalEstimatorDistance(featureI,featureJ) = optimalEstimatorDistance(featureI,featureJ) + ...
%                             featureDistance(optimalEstimatorMembers{featureI,featureJ}(featureK) , optimalEstimatorMembers{featureI,featureJ}(featureK+1));

                    optimalEstimatorDistance(featureI,featureJ) = nansum([optimalEstimatorDistance(featureI,featureJ) ...
                        featureDistance(optimalEstimatorMembers{featureI,featureJ}(featureK) , optimalEstimatorMembers{featureI,featureJ}(featureK+1)) ]);
                end

                displayCells{featureI,featureJ} = ['Feature ' num2str(featureI) '-' num2str(featureJ) ': '...
                    num2str(optimalEstimatorDistance(featureI,featureJ)) ' +- ' num2str(sqrt(optimalEstimatorVariance(featureI,featureJ)))];

            end
        end

                %%%%%%
%             for i =1:numFeatures-1
%                 for j = i+1:min(numFeatures,i+5)
%                     if optimalEstimatorMembers{i,j}(1) ~= i
%                         disp(num2str(optimalEstimatorMembers{i,j}));
%                     end
%                     if optimalEstimatorMembers{i,j}(end) ~= j
%                         disp(num2str(optimalEstimatorMembers{i,j}));
%                     end
%                 end
%             end
                %%%%%%

        optimalEstimatorStds = sqrt(optimalEstimatorVariance);

        mergableFeaturesEstimators = nan(numFeatures);
        mergableFeaturesEstimators(abs(optimalEstimatorDistance) < optimalEstimatorStds * confidenceInterval) = ...
            abs(optimalEstimatorDistance(abs(optimalEstimatorDistance) < optimalEstimatorStds * confidenceInterval));

        if ~any(mergableFeaturesEstimators(:))
            break;
        end

        while true


%             [~,featureMergeInd] = min(abs(mergableFeaturesEstimators(:)).*optimalEstimatorStds(:));
            [minVal,featureMergeInd] = min(mergableFeaturesEstimators(:));
            if isnan(minVal)
                break;
            end
            [mergeFeatureA,mergeFeatureB] = ind2sub([numFeatures numFeatures],featureMergeInd);

%                 estimatorTime(kymoNum) = toc(estimatorTimer)

%                 if mergeFeatureB < mergeFeatureA
%                     tempFeature = mergeFeatureB;
%                     mergeFeatureB = mergeFeatureA;
%                     mergeFeatureA = tempFeature;
%                 end
%                 numFeatures

            newFeature = merge_features(featuresCellArray_ordered,optimalEstimatorMembers{mergeFeatureA,mergeFeatureB},minOverlap);

            if newFeature == 0
                mergableFeaturesEstimators(mergeFeatureA,mergeFeatureB) = NaN;
                continue;
            end

%             featuresCellArray_merged = featuresCellArray_ordered;
            featuresCellArray_ordered{mergeFeatureA} = newFeature;
            featuresCellArray_ordered(mergeFeatureB) = [];


            featureDistance(:,mergeFeatureB) = [];
            featureDistance(mergeFeatureB,:) = [];
            featureDistanceVariance(:,mergeFeatureB) = [];
            featureDistanceVariance(mergeFeatureB,:) = [];
            featureOverlaps(:,mergeFeatureB) = [];
            featureOverlaps(mergeFeatureB,:) = [];

            [tempDistances, tempVariances, tempOverlaps] = calculate_feature_distances(featuresCellArray_ordered,mergeFeatureA,[],5);

            newDataMask = ~isnan(tempDistances);
            featureDistance(newDataMask) = tempDistances(newDataMask);
            featureDistanceVariance(newDataMask) = tempVariances(newDataMask);
            featureOverlaps(newDataMask) = tempOverlaps(newDataMask);

            for featureI = 1:numFeatures-1


                if featureI ~= mergeFeatureA

%                         for i = 1:numFeatures-2
%                             for j = i+1:numFeatures-1
%                                 if ~isempty(optimalEstimatorMembers{i,j})
%                                     if optimalEstimatorMembers{i,j}(1) ~= i
%                                         disp(num2str(optimalEstimatorMembers{i,j}));
%                                     end
%                                     if optimalEstimatorMembers{i,j}(end) ~= j
%                                         disp(num2str(optimalEstimatorMembers{i,j}));
%                                     end
%                                 end
%                             end
%                         end

                    minAI = min(featureI,mergeFeatureA);
                    maxAI = max(featureI,mergeFeatureA);

                    minBI = min(featureI,mergeFeatureB);
                    maxBI = max(featureI,mergeFeatureB);

                    numNewEstimatorMembersA = length(optimalEstimatorMembers{minAI,maxAI})-2; 
                    numNewEstimatorMembersB = length(optimalEstimatorMembers{minBI,maxBI})-2;

    %                     if numNewEstimatorMembersA > 0 && numNewEstimatorMembersB > 0
                    if numNewEstimatorMembersB > 0

%                             for i = 1:numFeatures-2
%                                 for j = i+1:numFeatures-1
%                                     if ~isempty(optimalEstimatorMembers{i,j})
%                                         if optimalEstimatorMembers{i,j}(1) ~= i
%                                             disp(num2str(optimalEstimatorMembers{i,j}));
%                                         end
%                                         if optimalEstimatorMembers{i,j}(end) ~= j
%                                             disp(num2str(optimalEstimatorMembers{i,j}));
%                                         end
%                                     end
%                                 end
%                             end

                        if numNewEstimatorMembersA > 0
                            optimalEstimatorMembers{minAI,maxAI}(end+numNewEstimatorMembersB) = optimalEstimatorMembers{minAI,maxAI}(end);
                            optimalEstimatorMembers{minAI,maxAI}(end-numNewEstimatorMembersB:end-1) = optimalEstimatorMembers{minBI,maxBI}(2:end-1);

%                                 for i = 1:numFeatures-2
%                                         for j = i+1:numFeatures-1
%                                             if ~isempty(optimalEstimatorMembers{i,j})
%                                                 if optimalEstimatorMembers{i,j}(1) ~= i
%                                                     disp(num2str(optimalEstimatorMembers{i,j}));
%                                                 end
%                                                 if optimalEstimatorMembers{i,j}(end) ~= j
%                                                     disp(num2str(optimalEstimatorMembers{i,j}));
%                                                 end
%                                             end
%                                         end
%                                     end

                        else if numNewEstimatorMembersA == 0

                               optimalEstimatorMembers{minAI,maxAI}(end+numNewEstimatorMembersB) = optimalEstimatorMembers{minAI,maxAI}(end);
                               optimalEstimatorMembers{minAI,maxAI}(end-numNewEstimatorMembersB:end-1) = optimalEstimatorMembers{minBI,maxBI}(2:end-1);

%                                    for i = 1:numFeatures-2
%                                         for j = i+1:numFeatures-1
%                                             if ~isempty(optimalEstimatorMembers{i,j})
%                                                 if optimalEstimatorMembers{i,j}(1) ~= i
%                                                     disp(num2str(optimalEstimatorMembers{i,j}));
%                                                 end
%                                                 if optimalEstimatorMembers{i,j}(end) ~= j
%                                                     disp(num2str(optimalEstimatorMembers{i,j}));
%                                                 end
%                                             end
%                                         end
%                                     end

                               else if numNewEstimatorMembersA < 0

                                       optimalEstimatorMembers{minAI,maxAI}(1) = minAI;
                                       optimalEstimatorMembers{minAI,maxAI}(2:numNewEstimatorMembersB+1) = optimalEstimatorMembers{minBI,maxBI}(2:end-1);
                                       optimalEstimatorMembers{minAI,maxAI}(end+1) = maxAI;

%                                            for i = 1:numFeatures-2
%                                                 for j = i+1:numFeatures-1
%                                                     if ~isempty(optimalEstimatorMembers{i,j})
%                                                         if optimalEstimatorMembers{i,j}(1) ~= i
%                                                             disp(num2str(optimalEstimatorMembers{i,j}));
%                                                         end
%                                                         if optimalEstimatorMembers{i,j}(end) ~= j
%                                                             disp(num2str(optimalEstimatorMembers{i,j}));
%                                                         end
%                                                     end
%                                                 end
%                                            end

                                   end
                            end

                        end

                       removalIdx = optimalEstimatorMembers{minAI,maxAI} == optimalEstimatorMembers{minAI,maxAI}(1);
                       removalIdx = removalIdx + (optimalEstimatorMembers{minAI,maxAI} == optimalEstimatorMembers{minAI,maxAI}(end));
                       removalIdx(1) = 0;
                       removalIdx(end) = 0;
                       optimalEstimatorMembers{minAI,maxAI}(logical(removalIdx)) = [];

                       optimalEstimatorMembers{minAI,maxAI} = unique(optimalEstimatorMembers{minAI,maxAI},'stable');

%                            if optimalEstimatorMembers{minAI,maxAI} ~= maxAI
%                                switchInd = find(optimalEstimatorMembers{minAI,maxAI} == maxAI);
%                                optimalEstimatorMembers{minAI,maxAI}([switchInd end]) = optimalEstimatorMembers{minAI,maxAI}([end switchInd]);
%                            end

%                            for i = 1:numFeatures-2
%                                 for j = i+1:numFeatures-1
%                                     if ~isempty(optimalEstimatorMembers{i,j})
%                                         if optimalEstimatorMembers{i,j}(1) ~= i
%                                             disp(num2str(optimalEstimatorMembers{i,j}));
%                                         end
%                                         if optimalEstimatorMembers{i,j}(end) ~= j
%                                             disp(num2str(optimalEstimatorMembers{i,j}));
%                                         end
%                                     end
%                                 end
%                             end

                    end

%                         for i = 1:numFeatures-2
%                             for j = i+1:numFeatures-1
%                                 if ~isempty(optimalEstimatorMembers{i,j})
%                                     if optimalEstimatorMembers{i,j}(1) ~= i
%                                         disp(num2str(optimalEstimatorMembers{i,j}));
%                                     end
%                                     if optimalEstimatorMembers{i,j}(end) ~= j
%                                         disp(num2str(optimalEstimatorMembers{i,j}));
%                                     end
%                                 end
%                             end
%                         end

                end

%                     for i = 1:numFeatures-2
%                         for j = i+1:numFeatures-1
%                             if ~isempty(optimalEstimatorMembers{i,j})
%                                 if optimalEstimatorMembers{i,j}(1) ~= i
%                                     disp(num2str(optimalEstimatorMembers{i,j}));
%                                 end
%                                 if optimalEstimatorMembers{i,j}(end) ~= j
%                                     disp(num2str(optimalEstimatorMembers{i,j}));
%                                 end
%                             end
%                         end
%                     end

            end

%                 for i = 1:numFeatures-2
%                     for j = i+1:numFeatures-1
%                         if ~isempty(optimalEstimatorMembers{i,j})
%                             if optimalEstimatorMembers{i,j}(1) ~= i
%                                 disp(num2str(optimalEstimatorMembers{i,j}));
%                             end
%                             if optimalEstimatorMembers{i,j}(end) ~= j
%                                 disp(num2str(optimalEstimatorMembers{i,j}));
%                             end
%                         end
%                     end
%                 end


            for featureI = 1:numFeatures-1
%                     for featureJ = 1:numFeatures-1
                for featureJ = featureI+1:numFeatures

                    membersMergeFeatureB = optimalEstimatorMembers{featureI,featureJ} == mergeFeatureB;

%                         if length(find(membersMergeFeatureB)) > 1
%                             a = 1;
%                         end
                    membersBeyondFeatureB = optimalEstimatorMembers{featureI,featureJ} > mergeFeatureB;
%                         membersBeyondFeatureB(membersMergeFeatureB) = 0;

                    optimalEstimatorMembers{featureI,featureJ}(membersMergeFeatureB) = mergeFeatureA;

%                         membersBeyondFeatureA = optimalEstimatorMembers{featureI,featureJ} > mergeFeatureA;
%                         optimalEstimatorMembers{featureI,featureJ}(membersBeyondFeatureA) = ...
%                             optimalEstimatorMembers{featureI,featureJ}(membersBeyondFeatureA)-1;

                    optimalEstimatorMembers{featureI,featureJ}(membersBeyondFeatureB) = ...
                        optimalEstimatorMembers{featureI,featureJ}(membersBeyondFeatureB)-1;

                end
            end

            optimalEstimatorMembers(mergeFeatureB,:) = [];
            optimalEstimatorMembers(:,mergeFeatureB) = [];

%                 for i = 1:numFeatures-2
%                     for j = i+1:numFeatures-1
%                         if ~isempty(optimalEstimatorMembers{i,j})
%                             if optimalEstimatorMembers{i,j}(1) ~= i
%                                 disp(num2str(optimalEstimatorMembers{i,j}));
%                             end
%                             if optimalEstimatorMembers{i,j}(end) ~= j
%                                 disp(num2str(optimalEstimatorMembers{i,j}));
%                             end
%                         end
%                     end
%                 end

%                 [a,replaceIdxs] = min(abs(mergableFeaturesEstimators(mergeFeatureA,mergeFeatureA+1:end)),abs(mergableFeaturesEstimators(mergeFeatureB,mergeFeatureA+1:end)));
%                 [a,replaceIdxs] = min(abs(squeeze(mergableFeaturesEstimators(mergeFeatureA,:))),squeeze(abs(mergableFeaturesEstimators(mergeFeatureB,:))));
%                 [b,replaceIdxs] = min(abs(mergableFeaturesEstimators(:,mergeFeatureA)),abs(mergableFeaturesEstimators(:,mergeFeatureB)));

            mergableFeaturesEstimators(mergeFeatureA,:) = min(mergableFeaturesEstimators(mergeFeatureA,:),mergableFeaturesEstimators(mergeFeatureB,:));
            mergableFeaturesEstimators(:,mergeFeatureA) = min(mergableFeaturesEstimators(:,mergeFeatureA),mergableFeaturesEstimators(:,mergeFeatureB));

            mergableFeaturesEstimators = min(mergableFeaturesEstimators,  mergableFeaturesEstimators.');
%                 mergableFeaturesEstimators(mergeFeatureA,mergeFeatureA+1:end) = min(mergableFeaturesEstimators(mergeFeatureA,:),mergableFeaturesEstimators(mergeFeatureB,:));
%                 mergableFeaturesEstimators(:,mergeFeatureA) = min(mergableFeaturesEstimators(:,mergeFeatureA),mergableFeaturesEstimators(:,mergeFeatureB));

            mergableFeaturesEstimators(logical(tril(ones(numFeatures,numFeatures),-1))) = NaN;
            mergableFeaturesEstimators(mergeFeatureA,mergeFeatureA) = NaN;
%                 mergableFeaturesEstimators(logical(tril(ones(numFeatures,numFeatures),0))) = NaN;

            mergableFeaturesEstimators(mergeFeatureB,:) = [];
            mergableFeaturesEstimators(:,mergeFeatureB) = [];

%                 mergableFeaturesEstimators(logical(eye(size(mergableFeaturesEstimators)))) = NaN;
%                 mergableFeaturesEstimators(mergeFeatureA,mergeFeatureA) = NaN;

%                 colorImg = get_color_map(featuresCellArray_ordered,size(kymo));
%                 figure, imagesc(colorImg);

            numFeatures = numFeatures-1;

        end

    end

%         for featureI = 1:numFeatures-1
%             for featureJ = featureI+1:min(numFeatures,featureI+5)
%                 disp(displayCells{featureI,featureJ});
%             end
%         end

%         toc;

%         kymoTime(kymoNum) = toc(kymoTimer)

    featureDistances = diag(optimalEstimatorDistance,1);
    [featureDistances, sorting] = order_features_by_dists(featureDistances);

%         featureDistanceVariances{kymoNum} = diag(featureDistanceVariance,1);
%         featureDistanceVariances{kymoNum} = featureDistanceVariances{kymoNum}(sorting);
    featureDistanceVariance = featureDistanceVariance(sorting,sorting);


    for featureI = 1:numFeatures-1
        [featureDistanceVariances(featureI), ~, ~] = graphshortestpath(sparse(featureDistanceVariance),featureI,featureI+1);
    end
%         featureDistanceVariances{kymoNum} = featureDistanceVariances{kymoNum}(sorting);

    featuresCellArray_processed = featuresCellArray_ordered(sorting);

%         kymoTime = toc(kymoTimer)
    
%     end
    
    fullTime = toc(fullTimer)
    
    

end
