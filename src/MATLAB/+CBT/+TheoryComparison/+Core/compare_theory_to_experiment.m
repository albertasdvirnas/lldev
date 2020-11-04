function comparisonResults = compare_theory_to_experiment(theoryStruct, experimentStruct, constantSettingsStruct, cacheResultsSubfolderPath, tryLoadPrecomputedResult, saveToPathTF,bs)
    if nargin < 4
        cacheResultsSubfolderPath = '';
    end
    if nargin < 5
        tryLoadPrecomputedResult = not(isempty(cacheResultsSubfolderPath));
    end
    if nargin < 6
        saveToPathTF =  not(isempty(cacheResultsSubfolderPath));
    end
    
    if nargin < 6
        bs = 1;
    end
    if tryLoadPrecomputedResult
        import CBT.TheoryComparison.Import.load_result_from_path;
        [found, comparisonResults] = load_result_from_path(cacheResultsSubfolderPath, theoryStruct, experimentStruct);
        if found
            return;
        end
    end
    
    import Fancy.Utils.extract_fields;
    [...
        stretchStrategy,...
        noOfStretchings,...
        optInterval,...
        lengthConstraint,...
        deltaCut,...
        psfSigmaWidth_nm,...
        psfSigmaWidth_bp,...
        pixelWidth_nm,...
        meanBpExt_pixels...
        ] = extract_fields(constantSettingsStruct, {...
            'stretchStrat',...
            'noOfStretchings',...
            'optInterval',...
            'lengthConstraint',...
            'deltaCut',...
            'psfWidth_nm',...
            'psfWidth_bp',...
            'nmPerPixel',...
            'pixelsPerBp'...
            });

    import CBT.TheoryComparison.get_struct_theory_curve_bpRes;
    theoryCurve_bpRes = get_struct_theory_curve_bpRes(theoryStruct,bs);
    
    import CBT.TheoryComparison.get_struct_experiment_curve_pxRes;
    experimentCurve_pxRes = get_struct_experiment_curve_pxRes(experimentStruct);
    experimentCurve_pxRes = zscore(experimentCurve_pxRes);
    
    import CBT.TheoryComparison.Import.get_experiment_bitmask;
    experimentCurveBitmask = get_experiment_bitmask(experimentStruct, deltaCut, psfSigmaWidth_nm, pixelWidth_nm); % get bitmask as if experiment
    if isempty(experimentCurveBitmask)
        warning('Empty curve');
    elseif not(any(experimentCurveBitmask))
        warning('No curve values were included');
    end
        
    import CBT.TheoryComparison.get_stretch_factors;
    stretchFactors = get_stretch_factors(length(theoryCurve_bpRes), length(experimentCurve_pxRes), meanBpExt_pixels, stretchStrategy, noOfStretchings, optInterval, lengthConstraint);
    
    import CBT.TheoryComparison.Core.compare_at_stretch_factors;
    [bestCC, meanCC, stdCC, bestStretchFactor]...
         = compare_at_stretch_factors(stretchFactors, theoryCurve_bpRes, psfSigmaWidth_bp, meanBpExt_pixels, experimentCurve_pxRes, experimentCurveBitmask);

    experimentStruct = rmfield(experimentStruct, 'experimentCurve_pxRes');
    comparisonResults = struct(...
        'bestCC', bestCC,...
        'meanCC', meanCC,...
        'stdCC', stdCC,...
        'structA', theoryStruct,...
        'structB', experimentStruct,...
        'bestStretchFactor', bestStretchFactor,...
        'stretchFactors', stretchFactors...
        );
    if saveToPathTF
        import CBT.TheoryComparison.Export.save_result_to_path;
        save_result_to_path(cacheResultsSubfolderPath, comparisonResults)
    end
end