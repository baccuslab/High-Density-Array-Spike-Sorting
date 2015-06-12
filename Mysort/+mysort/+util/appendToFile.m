
function appendToFile(filename, str)
  fh = fopen(filename, 'a+');
  fprintf(fh,str);
  fclose(fh);
end