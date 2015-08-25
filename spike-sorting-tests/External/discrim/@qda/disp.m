function disp(f)
%QDA/DISP Display QDA object.

%   Copyright (c) 1999 Michael Kiefte.

%   $Log$

ldet = f.ldet;
fprintf('    means: [%dx%d double]\n', size(f.means));
fprintf('    scale: [%dx%dx%d double]\n', size(f.scale));
fprintf('     ldet: [');
if length(f.ldet) < 10
  fprintf('%g ', ldet(1:end-1));
  fprintf('%g]\n', ldet(end));
else
  fprintf('1x%d double]\n', length(ldet));
end
fprintf('      est: ')
switch f.est
 case 0
  fprintf('0 (Unbiased estimator)\n');
 case 1
  fprintf('1 (Maximum-likelihood estimator)\n');
 otherwise
  fprintf('''t'' (t-parameter estimator)\n');
end
if f.est == 't'
  fprintf('       nu: %d\n', f.nu);
end
disp(f.classifier)




