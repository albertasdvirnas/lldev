function [] = plot_free_energy_contig_placements(...
        hAxis, ...
        refBarcode, ...
        placedContigBarcodes, ...
        placedContigMask, ...
        meanBpExt_pixels, ...
        branchIdx, ...
        sValsByBranch, ...
        numTotalOverlapsByBranch, ...
        contigPlacementOptionIdxsByBranch, ...
        shouldExportTxtTF, ...
        dataSampleName ...
        )

    sVal = sValsByBranch(branchIdx, 1);
    numTotalOverlap = numTotalOverlapsByBranch(branchIdx);
    contigPlacementOptionIdxs = contigPlacementOptionIdxsByBranch(branchIdx, :)';

    refBarcodeLabel = 'Experiment';

    placedContigLabels = cellfun(...
        @(contigIdx) ...
            sprintf('Contig %d', contigIdx), ...
            find(placedContigMask), ...
            'UniformOutput', false);

    refBarcodeLen_pixels = length(refBarcode);
    import CA.Core.gen_aligned_placement_mats;
    [~, refContigPlacementValsMat] = gen_aligned_placement_mats(refBarcodeLen_pixels, placedContigBarcodes, contigPlacementOptionIdxs);

    plotTitleStr = fprintf('Index: %d, P-value: %g, Overlap: %g', branchIdx, sVal, numTotalOverlap);
    import CA.UI.plot_contig_placements;
    kbpsPerPixel = 1/(meanBpExt_pixels * 1000);
    plot_contig_placements(...
        hAxis, ...
        plotTitleStr, ...
        refBarcode, ...
        refContigPlacementValsMat, ...
        kbpsPerPixel, ...
        refBarcodeLabel, ...
        placedContigLabels ...
    );

    % Save the plot to a text file
    if shouldExportTxtTF
        import CA.Export.export_contig_placements_txt;
        excludedContigsStr = '';
        excludedContigIdxs = find(not(placedContigMask));
        if not(isempty(excludedContigIdxs))
            excludedContigsStr = num2str(excludedContigIdxs);
        end
        export_contig_placements_txt(dataSampleName, refBarcode, refContigPlacementValsMat, meanBpExt_pixels, sVal, numTotalOverlap, excludedContigsStr);
    end
end