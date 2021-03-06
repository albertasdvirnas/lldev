function [] = add_dlc_menu(hMenuParent, tsDLT)
    hMenuCBC = uimenu( ...
        'Parent', hMenuParent, ...
        'Label', 'DLT');
    % add densely labeled consensus menu
    
    % run consensing
    import DLT.Consensus.run_consensing;
    uimenu( ...
        'Parent', hMenuCBC, ...
        'Label', 'Run Consensus Generation', ...
        'Callback', @(~, ~) run_consensing(tsDLT));
    
    % load consensus results
    import DLT.Consensus.Import.load_consensus_results;
    uimenu( ...
        'Parent', hMenuCBC, ...
        'Label', 'Load Consensus Results', ...
        'Callback', @(~, ~) load_consensus_results(tsDLT));
end