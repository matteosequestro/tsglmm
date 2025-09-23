

% Made this function because: (a) I never remember which one is a string
% and which one is a character; and (b) the default matlab functions have 
% stupid names why isn't it called char2str or str2char?
function x = changetext(x)

    if isa(x, 'char')
        x = convertCharsToStrings(x);
    elseif isa(x, 'string')
        x = convertStringsToChars(x);
    else
        warning("Conversion failed. Input must be either a character or a string")
    end


end

