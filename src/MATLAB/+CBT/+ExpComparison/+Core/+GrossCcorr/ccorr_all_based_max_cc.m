function [maxCC] = ccorr_all_based_max_cc(...
        barcodeA, ...
        barcodeB, ...
        longerBarcodeIsCircularTF, ...
        shouldRescaleTF ...
        )
    import CBT.ExpComparison.Core.GrossCcorr.ccorr_all;
    [betterDirCCs, ~, ~, ~] = ccorr_all(...
        barcodeA, ...
        barcodeB, ...
        longerBarcodeIsCircularTF, ...
        shouldRescaleTF ...
        );
    maxCC = max(betterDirCCs);
end