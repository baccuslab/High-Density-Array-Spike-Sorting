
function logToFile(filename, str)
  mysort.util.appendToFile(filename, [datestr(now) '  ' str '\n']);
end