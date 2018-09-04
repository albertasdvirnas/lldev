function bayerPattern = prompt_bayer_pattern()
    bayerPatternOptions = { ...
        'BGGR'; ...
        'GRBG'; ...
        'RGGB'; ...
        'GBRG'; ...
    };
    import Fancy.UI.FancyInput.dropdown_dialog;
    bayerPattern = lower(dropdown_dialog('Bayer Pattern Selection', 'Choose bayer pattern', bayerPatternOptions));
end