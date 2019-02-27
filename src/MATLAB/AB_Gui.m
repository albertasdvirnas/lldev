function output = AB_Gui(sets)
    % AB_Gui Autobarcoding main function

    %     Args:
    %         sets (struct): Input settings to the method
    % 
    %     Returns:
    %         output: output structure
    % 
    %     Example:
    
    output = [];
    
    % begin by loading default settings
    if nargin < 1
        [sets] = AB.Scripts.ab_sets();
        sets.moviesets.askformovie = 1;
    end
  
    % load all user selected settings
    sets = AB.UI.get_user_settings(sets);

    % these will save barcodes, backgrounds, and barcode display names
    barcodes = {};
    kymos = {};
    backgrounds = {};
    barcodeDisplayNames = {};
    
    % loop over movie file folder
    for idx = 1:length(sets.moviefilefold)
        
        % load one of the movies
        import AB.Processing.load_movie;
        movie = load_movie(sets.moviefilefold{idx},sets.filenames{idx});  
        
        % put it into a 3d matrix. Is that  necessary?
        movie3d = double(cat(3, movie{:}));
        
      %  slice = movie3d(:,:,1);
        
        % process movie by detecting the angle, rotating, detecting the
        % background mask, connected components, and use this to extract 
        import AB.Processing.preprocess_movie_for_kymo_extraction;
        [kymosMolEdgeIdxs, movieRot, rRot, cRot, ccIdxs, ccStructNoEdgeAdj, rotationAngle,bgvals] = preprocess_movie_for_kymo_extraction(movie3d, sets.preprocessing);

        % generate kymos
        import AB.Processing.generate_kymos_from_movie;
        abStruct = generate_kymos_from_movie(movie3d,kymosMolEdgeIdxs, rRot, cRot,sets.moviefilefold{idx},sets.kymo,sets.filenames{idx});
%         
%         % TODO: this is temporary hack to add noise to left and right. later change hca
%         % so that we don't need to detect edges there!
%         currRng = rng(); % so that current rng state can be restored
%         % temporarily set to produce predictable pseudorandom values for reproducibility
%         rng(rng(0, 'twister'));  
%         randVals = randn([size(abStruct.flattenedKymos{1},1), 50]) .*3* nanstd(movieRot(:)) + nanmean(movieRot(:));
%         rng(currRng); % restore rng state
%         for i=1:length(abStruct.flattenedKymos)
%             abStruct.flattenedKymos{i} = [randVals abStruct.flattenedKymos{i} randVals];
%         end
    
        kymos = [kymos; abStruct.flattenedKymos];
        barcodes = [barcodes; abStruct.barcodes];
        backgrounds = [backgrounds; abStruct.backgrounds];
        barcodeDisplayNames = [barcodeDisplayNames; abStruct.barcodeDisplayNames];
    end

    
    % here save kymo's to some temporary folder so we can filter out

    % generate consensus. here we can use kymos from all the movies
    import AB.Processing.compute_cluster_consensus;
    [lC, clusterMeanCenters, consensusInputs, consensusStructs] = compute_cluster_consensus(barcodes,backgrounds, barcodeDisplayNames, sets.consensus);
    
    
    import AB.UI.plot_consensus;
    hfig = plot_consensus(consensusStructs);
   
    %output  = movieProcessingResultsStruct;
end