function disp(f)
%SOFTMAX/DISP Display SOFTMAX object.

%   Copyright (c) 1999 Michael Kiefte.

%   $Id: disp.m,v 1.1 1999/06/04 18:50:50 michael Exp $
%   $Log: disp.m,v $
%   Revision 1.1  1999/06/04 18:50:50  michael
%   Initial revision
%

fprintf('  weights: [%dx%d sparse]\n', size(f.weights));
disp(f.classifier)