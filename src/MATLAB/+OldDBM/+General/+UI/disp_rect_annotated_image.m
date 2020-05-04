function [] = disp_rect_annotated_image(hAxis, imgArr, imgHeaderText, moleculeRectPositions)
    import OldDBM.General.UI.disp_img_with_header;
    disp_img_with_header(hAxis, imgArr, imgHeaderText);
    hold(hAxis, 'on');
    try
    cellfun(@(moleculeRectPosition) ...
        rectangle(...
            'Position', moleculeRectPosition, ...
            'LineWidth', 0.2, ...
            'EdgeColor', 'r'), ...
        moleculeRectPositions);
    catch
    end
    hold(hAxis, 'off');
end