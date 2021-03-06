function [theoryBarcodesZscaled_pxRes, theoryBarcodesUnscaled_pxRes, seqDescriptions, seqFilepaths, seqIdxsInFile, barcodeGenSettings] = quick_cbt_demo(barcodeGenSettings, fastaFilepaths, exportTsvs)
    % QUICK_CBT_DEMO
    %   Demo script for Competitive Binding Theory Curve generation
    %
    % Inputs:
    %   barcodeGenSettings
    %    (optional; CBT.ini settings prompt otherwise)
    %    struct containing fields for the following CBT/microscopy
    %    settings parameters:
    %      concNetropsin_molar
    %        concentration of Netropsin (molar)
    %      concYOYO1_molar
    %        concentration of YOYO-1 (molar)
    %      meanBpExt_nm
    %        length of basepair extension along channel (nm)
    %      psfSigmaWidth_nm
    %        width of sigma of point spread function (nm)
    %      pixelWidth_nm
    %        width of pixels of the CCD camera (nm)
    %   fastaFilepaths
    %     (optional; user is prompted by default if not provided)
    %      Cell array of the the filepaths for the fasta files to load
    %       and generate barcodes for
    %   exportTsvs
    %     (optional; false by default)
    %     whether user is prompted for a tsv location to save each curve
    %
    %  Authors:
    %    Saair Quaderi
    if (nargin < 1) || isempty(barcodeGenSettings)
        import CBT.get_default_barcode_gen_settings;
        barcodeGenSettings = get_default_barcode_gen_settings();

        import CBT.Import.Helpers.read_CBT_settings_struct;
    	cbtSettingsStruct = read_CBT_settings_struct();
        
        import Fancy.Utils.extract_fields;
        [ ...
            barcodeGenSettings.concNetropsin_molar, ...
            barcodeGenSettings.concYOYO1_molar ...
        ] = extract_fields(cbtSettingsStruct.cbtheory, {'NETROPSINconc', ...
                                                        'YOYO1conc'});
        [ ...
            barcodeGenSettings.meanBpExt_nm,...
            barcodeGenSettings.pixelWidth_nm,...
            barcodeGenSettings.psfSigmaWidth_nm,...
            barcodeGenSettings.stretchStrat,...
            barcodeGenSettings.noOfStretchings,...
            barcodeGenSettings.optInterval,...
            barcodeGenSettings.deltaCut,...
            barcodeGenSettings.lengthConstraint,...
            barcodeGenSettings.maxPairLengthDiffRelative,...
            barcodeGenSettings.maxPairLengthDiffAbsolute_nm...
            ] = extract_fields(...
            cbtSettingsStruct.converttheory, {...
            'nmPerBps',...
            'nmPerPixel',...
            'PSF_width',...
            'stretchStrat',...
            'noOfStretchings',...
            'optInterval',...
            'DeltaCut',...
            'lengthConstraint',...
            'maxPairLengthDiffRelative',...
            'maxPairLengthDiffAbsolute_nm'
            });                                         
        % warning('Using default barcode generation settings:');
    end
    
    if (nargin < 2) || isempty(fastaFilepaths)
        import NtSeq.Import.UI.try_prompt_nt_seq_filepaths;
        [~, fastaFilepaths] = try_prompt_nt_seq_filepaths([], true, false);
        if isempty(fastaFilepaths)
            error('No theory sequences were provided');
        end
    end
    
    if (nargin < 3) || isempty(exportTsvs)
        exportTsvs = false;
    end
    
    paramTableVarNames = matlab.lang.makeValidName({...
        'Netropsin concentration (molar)', ...
        'YOYO-1 concentration (molar)', ...
        'Basepair extension len (nm)', ...
        'PSF sigma width (nm)', ...
        'Pixel width (nm)'}');
    paramsTable = table( ...
        barcodeGenSettings.concNetropsin_molar, ...
        barcodeGenSettings.concYOYO1_molar, ...
        barcodeGenSettings.meanBpExt_nm, ...
        barcodeGenSettings.psfSigmaWidth_nm, ...
        barcodeGenSettings.pixelWidth_nm, ...
        'VariableNames', paramTableVarNames);
    disp(paramsTable);
    

    hFig = figure('Name', 'Competitive Binding Theory');
    hParent = uipanel('Parent', hFig);
    
    % -- import fasta data
    import NtSeq.Import.import_fasta_nt_seqs;
    [ntSeqs, seqFastaHeaders, seqFilepaths, seqIdxsInFile] = import_fasta_nt_seqs(fastaFilepaths);
    
    % compute competitive binding theory barcode
    import CBT.Core.gen_unscaled_cbt_barcodes;
    theoryBarcodesUnscaled_pxRes = gen_unscaled_cbt_barcodes(ntSeqs, barcodeGenSettings);
    
    %   theoryBarcodes_pxRes_zscaled
    %      Cell array of the barcodes after point spread function
    %       convolution and sampling to pixel resolution and
    %      rescaling of mean and variance to zero and one respectively
    
    % compute rescaled barcodes with a mean of 0 and variance of 1
    theoryBarcodesZscaled_pxRes = cellfun(@zscore, theoryBarcodesUnscaled_pxRes, 'UniformOutput', false);

    
    
    import CBT.UI.plot_curves_in_tabs;
    
    [~, seqFilenamesSansExt] = cellfun(@fileparts, ...
        seqFilepaths, ...
        'UniformOutput', false);
    tabTitles = seqFilenamesSansExt;
    seqDescriptions = cellfun( ...
        @(fastaFilenameSansExt, fastaHeader) ...
            sprintf('Barcode for %s\n%s', ...
                fastaFilenameSansExt, ...
                strtrim(fastaHeader(max([0, strfind(fastaHeader, '|')]) + 1:end))), ...
        seqFilenamesSansExt, seqFastaHeaders, ...
        'UniformOutput', false);
    hTabGroup = uitabgroup('Parent', hParent);
    plot_curves_in_tabs(hTabGroup, theoryBarcodesUnscaled_pxRes, tabTitles, seqDescriptions);
    
    
    if exportTsvs
        import CBT.UI.prompt_tsv_output_filepaths;
        import CBT.Export.export_theory_curve_tsvs;

        [tsvFilepaths] = prompt_tsv_output_filepaths(seqFilepaths, seqFilenamesSansExt);
        export_theory_curve_tsvs(theoryBarcodesZscaled_pxRes, seqFilenamesSansExt, tsvFilepaths);
    end
end