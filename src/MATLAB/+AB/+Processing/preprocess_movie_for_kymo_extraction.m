function [kymoMolEdgeIdxs, movieRot, rRot, cRot, ccIdxs, ccStructNoEdgeAdj, rotationAngle,bgvals] = preprocess_movie_for_kymo_extraction(movieIn, preprocessSettings)
    % preprocess_movie_for_kymo_extraction
    
    % :param movieIn: input movie, should be 3d array
    % :param preprocessSettings: setting to process the movie
    %
    % :returns: kymoMolEdgeIdxs, movieRot, rRot, cRot, ccIdxs, ccStructNoEdgeAdj, rotationAngle
    
    % rewritten by Albertas Dvirnas
    
	%import AB.Processing.get_rotation_angle;
    import OldDBM.MoleculeDetection.get_angle
    
    import AB.Core.rotate_movie;
    import AB.Processing.get_foreground_mask_movie;
    import AB.Core.get_filtered_cc_struct;

    %% Registration
    
    % might need to do registration for some of the movies, add this as an
    % option, perhaps this is a good place for that?
    
    %% Derect angle
    
    %[rotationAngle] = get_rotation_angle(movieIn, preprocessSettings.rotation);
    rotationAngle = get_angle(movieIn);

    %% Rotate movie
    [movieRot, rRot, cRot, segmentFrameRot] = rotate_movie(movieIn, -rotationAngle);
    
     % add a consistency check to see if the movie rotated correctly.
	I1 = double(nanmean(movieIn, 3));
    I2 = double(nanmean(movieRot, 3));
    if max(nanmean(I1)) > max(nanmean(I2))
       warning('The angle might not have been computed correctly, since the highest peak is smaller after rotation'); 
    end

    %% Assign nan's to all the pixels nonrepresented in the movie
    movieRot(repmat(~segmentFrameRot, [1, 1, size(movieRot, 3)])) = NaN;
    
    %% compute a foreground mask
    [fgMaskMov,bgvals] = get_foreground_mask_movie(movieRot, preprocessSettings.foregroundMasking);
   
%     figure,ax1=subplot(1,2,1);imshow(fgMaskMov(:,:,5),[]);
%     ax2=subplot(1,2,2);imshow(movieRot(:,:,5),[]);
%     linkaxes([ax1 ax2]);
    
    % ------
    fprintf('Conducting filtered connected component analysis on molecules\n');
    tic
    
    % filter out components that aren't present in at least this many
    % frames. They should be present in all the frames, otherwise there is too much
    % movement
    minFramePresence = size(fgMaskMov,3); 
    % minimum length in x direction
    minLen = 30;
    [ccStructNoEdgeAdj] = get_filtered_cc_struct(fgMaskMov, minFramePresence, minLen, segmentFrameRot);
    
    toc
    fprintf('Conducted filtered connected component analysis on molecules\n');

    % ------
    fprintf('Detecting kymograph edges\n');
    tic

    import AB.Core.detect_kymo_mol_edge_idxs;
    [kymoMolEdgeIdxs, ccIdxs] = detect_kymo_mol_edge_idxs(ccStructNoEdgeAdj, preprocessSettings.kymoEdgeDetection);

    toc
    fprintf('Detected kymograph edges\n');

end
