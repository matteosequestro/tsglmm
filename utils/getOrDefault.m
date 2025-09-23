function val = getOrDefault(s, field, defaultVal)
    if isfield(s, field)
        val = s.(field);
    else
        val = defaultVal;
    end
end