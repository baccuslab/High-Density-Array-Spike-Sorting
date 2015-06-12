function display(f)
%LDA/DISPLAY Display an LDA object.

%   Copyright (c) 1999 Michael Kiefte.

%   $Log$

fprintf('\n%s =\n\n', inputname(1));
fprintf('\tLinear discriminant analysis\n');
switch f.est
 case 0
  fprintf('\tMaximum likelihood estimator\n');
 case 1
  fprintf('\tBias corrected estimator\n');
 otherwise
  fprintf('\tStudent-t parameter estimator on %d degrees of freedom\n', f.nu);
end
fprintf('\t%d variates; %d classes\n\n', size(f.means'))
disp(f)
fprintf('\n');
