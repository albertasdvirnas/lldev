function comAlignedKymo = center_of_mass_align_kymo(rawKymo, centerOfMass)
    [numRows, numCols] = size(rawKymo);
    meanCenterOfMass = centerOfMass(1);%mean(centerOfMass);
    comAlignedKymo = zeros(numRows, numCols);
    for rowNum = 1:numRows
        singleKymoFrame = rawKymo(rowNum, :);
        singleFrameKymoNew = circshift(singleKymoFrame', [-round(centerOfMass(rowNum) - meanCenterOfMass) 0]);
        comAlignedKymo(rowNum, :) = singleFrameKymoNew;
    end
end