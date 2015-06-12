function display(f)
%SOFTMAX/DISPLAY Display SOFTMAX object.

%   Copyright (c) 1999 Michael Kiefte.

%   $Id: display.m,v 1.1 1999/06/04 18:50:50 michael Exp $
%   $Log: display.m,v $
%   Revision 1.1  1999/06/04 18:50:50  michael
%   Initial revision
%

weights = f.weights;
p = full(sum(~any(weights)))-1;
g = diff(size(weights))+1;
fprintf('\n%s =\n\n', inputname(1));
fprintf('\tMultinomial feed-forward neural network\n');
fprintf('\t%d variates; %d classes\n\n', p, g);
disp(f)
fprintf('\n');