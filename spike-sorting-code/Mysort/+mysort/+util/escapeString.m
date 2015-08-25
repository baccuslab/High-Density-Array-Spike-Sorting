
function str = escapeString(str)
    str = regexprep(str, '\\', '\\\\');
    