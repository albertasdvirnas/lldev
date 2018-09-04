function [hcaSessionStruct] = get_theory_to_exp_p_values_long(hcaSessionStruct, sets )
    % the main function to compare fragments of human chromosome vs. theory
    
    % input hcaSessionStruct, sets 
    
    % output comparisonStructure? (move it outside hcaSessionStruct?)
        
    disp('Starting computing exp to theory p-values...')
    tic

    % barcodeGenSettings = hcaSessionStruct.theoryGen.sets;
  %  import CA.CombAuc.Core.Zeromodel.generate_random_sequences;
    %import CA.CombAuc.Core.Comparison.compute_correlation;
    import CA.CombAuc.Core.Comparison.generate_evd_par;
   % import CA.CombAuc.Core.Comparison.cc_fft;

    % double check if this is the one needed
    if  sets.skipNullModelChoice == 0
        sets.nullModelPath = uigetdir(pwd,'Select folder with the pre-computed null model');
        addpath(genpath( sets.nullModelPath));    
    end
         
    meanFFTest = load(strcat([sets.nullModelPath '/meanFFT.mat']));
    meanFFTest = meanFFTest.meanFFTEst;   
    
    
    
    meanBpExt_nm =sets.barcodeGenSettings.meanBpExt_nm; % from theory lengths
    psfSigmaWidth_bps = sets.barcodeConsensusSettings.psfSigmaWidth_nm/meanBpExt_nm;
    sets.untrustedBp = round(sets.barcodeConsensusSettings.deltaCut*psfSigmaWidth_bps);
    pxPerBp = meanBpExt_nm/sets.barcodeConsensusSettings.prestretchPixelWidth_nm;
    bpPerPx = sets.barcodeConsensusSettings.prestretchPixelWidth_nm/meanBpExt_nm;

    refBarcode= hcaSessionStruct.theoryGen.theoryBarcodes{1};
    refBitmask= hcaSessionStruct.theoryGen.bitmask{1};
    sets.askForPvalueSettings=1;
    % don't hardcode this
    if sets.askForPvalueSettings == 1
        sets.contigSettings.numRandBarcodes = 1000;
        sets.pvaluethresh = 0.01;
    end
    
    if sets.prestretchMethod == 0
        % if the barcodes were not pre-stretched, it means that we need to
        % consider each on of them separately
        lengths = [hcaSessionStruct.lengths mean(hcaSessionStruct.lengths)];
        
        coefs = cell2mat(cellfun(@(x) x.maxcoef(1),hcaSessionStruct.comparisonStructure,'UniformOutput',0));
        pValueMatrix= ones(1,length(coefs));
        pValueMatrixFiltered = ones(1,length(coefs));
        rsq = ones(1,length(coefs));
       % pValueMatrix = CA.CombAuc.Core.Comparison.compute_p_value(coefs,evdPar,'exact'); 

        stretchFactors = sets.barcodeConsensusSettings.stretchFactors;

        % make this in a nice function
        for ii=1:length(coefs)
             import CA.CombAuc.Core.Zeromodel.generate_random_sequences;
            [ randomSequences ] = generate_random_sequences(2*round(lengths(refBarcode)*bpPerPx),sets.contigSettings.numRandBarcodes,meanFFTest,psfSigmaWidth_bps ,'phase',pxPerBp);
            
            ccMax = ones(length(stretchFactors),sets.contigSettings.numRandBarcodes);

            for i=1:length(randomSequences)
                randomSeq = randomSequences{i}(1:round(length(randomSequences{i})/2));
                %randomBit = ones(1,length(randomSeq));
                for j=1:length(stretchFactors)
                    barC = interp1(randomSeq, linspace(1,length(randomSeq),length(randomSeq)*stretchFactors(j)));
                    [xcorrs,~,~] =  CA.CombAuc.Core.Comparison.get_cc_fft(barC,randomSeq,ones(1,length(barC)),ones(1,length(randomSeq)));
                    ccMax(j,i) = max(xcorrs(:));
                end  
            end
            hcaSessionStruct.pValueResults.ccMax = ccMax;
            if length(stretchFactors) > 1
                ccMax = max(ccMax);
            end

           % for this method we use standart EVD generation (see contig paper)
            import CA.CombAuc.Core.Comparison.generate_evd_par;
            [~,rsq(ii),evdPar] = generate_evd_par( ccMax,[],length(randomSequences{1})/2, 'exact2' );


            pValueMatrix(ii) = CA.CombAuc.Core.Comparison.compute_p_value(coefs(ii),evdPar,'exact'); 

            if sets.filterSettings.filter==1
                coefs2 = cell2mat(cellfun(@(x) x.maxcoef(1),hcaSessionStruct.comparisonStructureFiltered,'UniformOutput',0));
                pValueMatrixFiltered(ii) = CA.CombAuc.Core.Comparison.compute_p_value(coefs2(ii),evdPar,'exact'); 
            end
            
        end
        hcaSessionStruct.pValueResults.evdPar = evdPar;
        hcaSessionStruct.pValueResults.pValueMatrix = pValueMatrix;
        hcaSessionStruct.pValueResults.rsq = rsq;
        hcaSessionStruct.pValueResults.pValueMatrixFiltered = pValueMatrixFiltered;     
    else
        % otherwise we compure random barcode for all at the same time

        % compute short random barcodes
        import CA.CombAuc.Core.Zeromodel.generate_random_sequences;
        [ randomSequences ] = generate_random_sequences(2*round(mean(hcaSessionStruct.lengths)*bpPerPx-2*sets.untrustedBp),sets.contigSettings.numRandBarcodes,meanFFTest,psfSigmaWidth_bps ,'phase',pxPerBp);
    %    toc
        stretchFactors = sets.barcodeConsensusSettings.stretchFactors;
        ccMax = ones(length(stretchFactors),sets.contigSettings.numRandBarcodes);

      %  xx =[];
        for i=1:length(randomSequences)
           % i
            randomSeq = randomSequences{i}(1:round(length(randomSequences{i})/2));
            randomBit = ones(1,length(randomSeq));
            for j=1:length(stretchFactors)
                barC = interp1(randomSeq, linspace(1,length(randomSeq),length(randomSeq)*stretchFactors(j)));
                [xcorrs,~,~] =  CA.CombAuc.Core.Comparison.get_cc_fft(barC,refBarcode,ones(1,length(barC)),refBitmask);
                ccMax(j,i) = max(xcorrs(:));
            end  
          %  [xcorrs,~,~] =  CA.CombAuc.Core.Comparison.get_cc_fft(randomSeq,refBarcode,ones(1,length(randomSeq)),refBitmask);
           %  ccMax(i) = max(xcorrs(:));
          %  xx = [xx;xcorrs(:)];
        end

        hcaSessionStruct.pValueResults.ccMax = ccMax;

        if length(stretchFactors) > 1
            ccMax = max(ccMax);
        end

           % ccMax(ccMax>1) = 1-1E-5;
       %    tic
       
       % for this method we use standart EVD generation (see contig paper)
        import CA.CombAuc.Core.Comparison.generate_evd_par;
        [~,rsq,evdPar] = generate_evd_par( ccMax,[],length(randomSequences{1})/2, 'exact2' );

        coefs = cell2mat(cellfun(@(x) x.maxcoef(1),hcaSessionStruct.comparisonStructure,'UniformOutput',0));

        pValueMatrix = CA.CombAuc.Core.Comparison.compute_p_value(coefs,evdPar,'exact'); 

        hcaSessionStruct.pValueResults.evdPar = evdPar;
        hcaSessionStruct.pValueResults.pValueMatrix = pValueMatrix;
        hcaSessionStruct.pValueResults.rsq = rsq;

        if sets.filterSettings.filter==1
            coefs = cell2mat(cellfun(@(x) x.maxcoef(1),hcaSessionStruct.comparisonStructureFiltered,'UniformOutput',0));
            pValueMatrix = CA.CombAuc.Core.Comparison.compute_p_value(coefs,evdPar,'exact'); 
            hcaSessionStruct.pValueResults.pValueMatrixFiltered = pValueMatrix;     
        end
        toc
        disp('Finished computing exp to theory p-values (method 1)...')
    end
	% Comparison.compare_distribution_to_data( ccMax(:), evdPar, 'exactfull' )
end

