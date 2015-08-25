function display(f)
%LOGDA/DISPLAY Display LOGDA object.

%   Copyright (c) 1999 Michael Kiefte.

%   $Log$

beta = f.coefs;
g = size(beta, 1) + 1;
p = size(beta, 2) - 1;
fprintf('\n%s =\n\n', inputname(1));
fprintf('\tLogistic discriminant analysis\n');
fprintf('\t%d variates; %d classes\n\n', p, g);
disp(f)
fprintf('\n');
