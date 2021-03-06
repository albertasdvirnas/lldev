function [] = SVD_Gui()
    % SVD_GUI - Structural Variation Detection (SVD) GUI (New) 

    hFig = figure(...
        'Name', 'Structural Variation Detection', ...
        ...%'Units', 'normalized', ...
        'OuterPosition', [0.05 0.05 650 500], ...
        'NumberTitle', 'off', ...
        'MenuBar', 'none' ...
        );
    hPanelSVD = uipanel('Parent', hFig);

    import Fancy.UI.FancyTabs.TabbedScreen;
    tsSVD = TabbedScreen(hPanelSVD);
    
    import SVD.UI.run_svd;
    run_svd(tsSVD);
end