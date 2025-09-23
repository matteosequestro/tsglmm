function new_text = switchText(text)

if isa(text, 'string')
    new_text = convertStringsToChars(text);

elseif isa(text, "char")
    new_text = convertCharsToStrings(text);
end


end