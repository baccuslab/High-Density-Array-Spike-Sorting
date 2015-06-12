function display(f)
%CLASSIFIER/DISPLAY Display CLASSIFIER object.

%   Copyright (c) 1999 Michael Kiefte.

%   $Log$

p = size(f.range, 2);
g = length(f.counts);
fprintf('\n%s =\n\n', inputname(1));
fprintf('\tClassifier object\n');
fprintf('\t%d variates; %d classes\n\n', p, g);
disp(f)
fprintf('\n');
