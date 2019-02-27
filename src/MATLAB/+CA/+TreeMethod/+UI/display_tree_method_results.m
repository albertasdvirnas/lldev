function [] = display_tree_method_results(ts, scaledRefBarcode, croppedContigBarcodes, passesMinimalLenMask, bpsPerPixel, dataSampleName, numPixelsTrimmed, contigOrderingsMat, startMat, sValsByBranch, coverageByBranch, flippedMat, sValsHistMat)
    hTabResults = ts.create_tab('Tree Method Results');
    hPanelResults = uipanel('Parent', hTabResults);
    hAxis = axes(...
        'Parent', hPanelResults, ...
        'Units', 'normalized', ...
        'Position', [0.1 0.25 0.83 0.7]);


    % Plotting the histograms
    import CA.UI.plot_hists;
    plot_hists(hAxis, sValsHistMat);



    [~, minSvalIdx] = min(sValsHistMat(3,:));

    % Create the buttons
    hEditBranchIdxs = uicontrol('Parent', hPanelResults,...
        'Units','normalized',...
        'Position',[0.77 0.07 0.18 0.07],...
        'Style','edit',...
        'String',num2str(minSvalIdx));
    hShouldSaveTxtTFRadio = uicontrol('Parent', hPanelResults,...
        'Style','radiobutton',...
        'String','Save plot as text file',...
        'Value',0,...
        'Units','normalized',...
        'Position',[0.4 0.02 0.2 0.03]);
    import CA.Export.export_contig_assembly_tree_txt;
    uicontrol('Parent', hPanelResults,...
        'Units','normalized',...
        'Position',[0.29 0.07 0.18 0.07],...
        'String','Save txt file',...
        'Callback', @(~, ~) ...
            export_contig_assembly_tree_txt( ...
                sValsByBranch, ...
                coverageByBranch, ...
                contigOrderingsMat, ...
                startMat, ...
                cellfun(@str2double, strsplit(get(hEditBranchIdxs,'String'), ',')), ...
                croppedContigBarcodes, ...
                num2str(find(not(passesMinimalLenMask))), ... 
                numPixelsTrimmed, ...
                dataSampleName ...
            ));
    import CA.Export.export_contig_assembly_hists_fig;
    uicontrol('Parent', hPanelResults,...
        'Units','normalized',...
        'Position',[0.05 0.07 0.18 0.07],...
        'String','Save histogram',...
        'Callback',@(~, ~) export_contig_assembly_hists_fig(sValsHistMat));

    import CA.TreeMethod.UI.plot_tree_method_branch_results
    uicontrol('Parent', hPanelResults,...
        'Units', 'normalized',...
        'Position', [0.53 0.07 0.18 0.07],...
        'String','Plot branch',...
        'Callback', @ (~, ~) ...
            feval(@(branchIdxs) ...
                plot_tree_method_branch_results( ...
                    get_branches_tabscreen(hTabResults), ...
                    dataSampleName, ...
                    scaledRefBarcode, ...
                    croppedContigBarcodes, ...
                    (bpsPerPixel ./ 1000), ...
                    branchIdxs, ...
                    startMat, ...
                    flippedMat, ...
                    contigOrderingsMat, ...
                    get(hShouldSaveTxtTFRadio,'Value') ...
                ), cellfun(@str2double, strsplit(get(hEditBranchIdxs,'String'), ','))) ...
            );
    function tsBranches = get_branches_tabscreen(hTabResults)
        persistent localTsBranches;
        if not(isempty(localTsBranches))
            tsBranches = localTsBranches;
            return;
        end
        
        hTabBranches = hTabResults.create_tab('Tree Method Results');
        hPanelBranches = uipanel('Parent', hTabBranches);
        
        import Fancy.UI.FancyTabs.TabbedScreen;
        tsBranches = TabbedScreen(hPanelBranches);
        localTsBranches = tsBranches;
    end
end