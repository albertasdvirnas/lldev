function [] = get_parameters_ui(ts, continue_callback)
    import Fancy.UI.FancyPositioning.get_pixel_height;
    import Fancy.Utils.extract_fields;

    tabTitle = 'Theory Parameters';

    hTheoryParamsTab = ts.create_tab(tabTitle);
    hTheoryParamsPanel = uipanel(hTheoryParamsTab, 'Position', [0, 0, 1, 1]);
    ts.select_tab(hTheoryParamsTab);

    import CBT.Import.Helpers.read_CBT_settings_struct;
    cbtSettingsStruct = read_CBT_settings_struct();

    [ ...
        concNetropsin_molar,...
        concYOYO1_molar...
        ] = extract_fields(...
        cbtSettingsStruct.cbtheory, {...
        'NETROPSINconc',...
        'YOYO1conc'...
        });
    [ ...
        meanBpExt_nm,...
        pixelWidth_nm,...
        psfSigmaWidth_nm,...
        stretchStrat,...
        noOfStretchings,...
        optInterval,...
        deltaCut,...
        lengthConstraint,...
        maxPairLengthDiffRelative,...
        maxPairLengthDiffAbsolute_nm,...
        bindingSequence...
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
        'maxPairLengthDiffAbsolute_nm',...
        'bindingSequence'
        });
    %todo: default defaults
    checkCacheForIntensityCurves = true;
    saveToCache = true;

    % Column names and column format
    columnName = {'Sequences',...
        '<html>Binding<br/>sequence</html>',...
        '<html>Netropsin<br/>concentration</html>',...
        '<html>YOYO-1<br/>concentration</html>',...
        '<html>nm<br/>per<br/>basepair</html>',...
        '<html>pixel<br/>width<br/>(nm)</html>',...
        '<html>PSF sigma<br/>width<br/>(nm)</html>',...
        '<html>Stretch Strategy</html>',...
        '<html># of stretchings<br/>[strats 3/4]</html>',...
        '<html>Stretch interval<br/>(optInterval)<br/>[stats 3/4]</html>',...
        '<html>Delta cut<br/>[x PSF width]</html>',...
        '<html>Max Pair<br/>Length Diff<br/>Relative</html>',...
        '<html>Max Pair<br/>Length Diff<br/>Absolute (nm)</html>',...
        '<html>Length<br/>constraint<br/>[strats 3/4]</html>',...
        '<html>Check<br/>Cache</html>',...
        '<html>Save<br/>to<br/>Cache</html>'};
    stretchStrats = {'1: none', '2: equal length', '3: optimal length', '4: hybrid of 2 & 3'};
    columnFormat = {{'*'},...
        'char',...
        'numeric',...
        'numeric',...
        'numeric',...
        'numeric',...
        'numeric',...
        stretchStrats,...
        'numeric',...
        'numeric',...
        'numeric',...
        'numeric',...
        'numeric',...
        'logical',...
        'logical',...
        'logical'};
    columnEditable = [...
        false,...
        true,...
        true,...
        true,...
        true,...
        true,...
        true,...
        true,...
        true,...
        true,...
        true,...
        true,...
        true,...
        true,...
        false];
    columnUneditable = columnEditable & false;

    % Define the data
    inData = {'*',...
        bindingSequence,...
        concNetropsin_molar,...
        concYOYO1_molar,...
        meanBpExt_nm,...
        pixelWidth_nm,...
        psfSigmaWidth_nm,...
        stretchStrats{stretchStrat},...
        noOfStretchings,...
        optInterval,...
        deltaCut,...
        maxPairLengthDiffRelative,...
        maxPairLengthDiffAbsolute_nm,...
        logical(lengthConstraint),...
        checkCacheForIntensityCurves,...
        saveToCache};

    %TODO: validate default inputs

    % Create the uitable
    hParamTable = uitable(hTheoryParamsPanel,...
        'Data', inData,... 
        'ColumnName', columnName,...
        'ColumnFormat', columnFormat,...
        'ColumnEditable', columnEditable,...
        'RowName',[]);
    hButtonContinue = uicontrol(hTheoryParamsPanel, 'Style', 'pushbutton', 'String', 'Continue to Curve Generation');

    parentHeightPx = get_pixel_height(get(hParamTable, 'Parent'));
    % Set width and height
    PADDING_HEIGHT_PX = 20;
    tableWidthPx = hParamTable.Extent(3);
    tableHeightPx = hParamTable.Extent(4);
    tableBottomPosPx = parentHeightPx - tableHeightPx - PADDING_HEIGHT_PX - 50;
    hParamTable.Position(3) = tableWidthPx;
    hParamTable.Position(4) = tableHeightPx;
    hParamTable.Position(2) = tableBottomPosPx;

    buttonWidthPx = 200;
    buttonHeightPx = 30;
    hButtonContinue.Position(3) = buttonWidthPx;
    hButtonContinue.Position(4) = buttonHeightPx;
    hButtonContinue.Position(2) = tableBottomPosPx - buttonHeightPx - PADDING_HEIGHT_PX;


    function [...
            bindingSequence,...
            concNetropsin_molar,...
            concYOYO1_molar,...
            meanBpExt_nm,...
            pixelWidth_nm,...
            psfSigmaWidth_nm,...
            stretchStrat,...
            noOfStretchings,...
            optInterval,...
            deltaCut,...
            maxPairLengthDiffRelative,...
            maxPairLengthDiffAbsolute_nm,...
            lengthConstraint,...
            checkCacheForIntensityCurves,...
            saveToCache] = get_parameters()
        outData = get(hParamTable, 'Data');
        [...
            bindingSequence,...
            concNetropsin_molar,...
            concYOYO1_molar,...
            meanBpExt_nm,...
            pixelWidth_nm,...
            psfSigmaWidth_nm,...
            stretchStrat,...
            noOfStretchings,...
            optInterval,...
            deltaCut,...
            maxPairLengthDiffRelative,...
            maxPairLengthDiffAbsolute_nm,...
            lengthConstraint,...
            checkCacheForIntensityCurves,...
            saveToCache] = outData{2:end};
        stretchStrat = floor(str2double(stretchStrat(1)));
    end
%             
%             %TODO: validate input changes
%             function on_cell_edit(~, ~)
%                 % validate
%             end
%             iptaddcallback(paramTableHandle, 'CellEditCallback', {@on_cell_edit});
%             

    function continue_click(~, ~)
        [bindingSequence,...
            concNetropsin_molar,...
            concYOYO1_molar,...
            meanBpExt_nm,...
            pixelWidth_nm,...
            psfSigmaWidth_nm,...
            stretchStrat,...
            noOfStretchings,...
            optInterval,...
            deltaCut,...
            maxPairLengthDiffRelative,...
            maxPairLengthDiffAbsolute_nm,...
            lengthConstraint,...
            checkCacheForIntensityCurves,...
            saveToCache] = get_parameters();

        if stretchStrat < 3
            noOfStretchings = 1;
            optInterval = 0;
        end

        paramsStruct = struct(...
            'bindingSequence', bindingSequence,...
            'NETROPSINconc', concNetropsin_molar,...
            'YOYO1conc', concYOYO1_molar,...
            'nmPerBp', meanBpExt_nm,...
            'nmPerPixel', pixelWidth_nm,...
            'pixelsPerBp', (meanBpExt_nm/pixelWidth_nm),...
            'psfWidth_nm', psfSigmaWidth_nm,...
            'psfWidth_bp', psfSigmaWidth_nm/meanBpExt_nm,...
            'stretchStrat', stretchStrat,...
            'noOfStretchings', noOfStretchings,...
            'optInterval', optInterval,...
            'deltaCut', deltaCut,...
            'maxPairLengthDiffRelative', maxPairLengthDiffRelative,...
            'maxPairLengthDiffAbsolute_nm', maxPairLengthDiffAbsolute_nm,...
            'lengthConstraint', lengthConstraint,...
            'checkCacheForIntensityCurves', checkCacheForIntensityCurves,...
            'saveToCache', saveToCache...
            );

        continue_callback(paramsStruct);
        set(hParamTable, 'ColumnEditable', columnUneditable);
    end
    iptaddcallback(hButtonContinue, 'Callback', {@continue_click});
end