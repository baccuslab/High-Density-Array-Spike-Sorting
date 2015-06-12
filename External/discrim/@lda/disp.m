function disp(f)
%LDA/DISP Display LDA object.

%   Copyright (c) 1999 Michael Kiefte.

%   $Log$

fprintf('    means: [%dx%d double]\n', size(f.means));
fprintf('    scale: [%dx%d double]\n', size(f.scale));
fprintf('      est: ')
switch f.est
 case 1
  fprintf('1 (Maximum-likelihood estimator)\n');
 case 0
  fprintf('0 (Unbiased estimator)\n');
 otherwise
  fprintf('''t'' (t-parameter estimator)\n');
end
if f.est == 't'
  fprintf('       nu: %d\n', f.nu);
end
disp(f.classifier)






