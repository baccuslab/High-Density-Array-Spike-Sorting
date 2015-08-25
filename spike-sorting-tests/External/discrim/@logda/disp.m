function disp(f)
%LOGDA/DISP Display LOGDA object.

%   Copyright (c) 1999 Michael Kiefte.

%   $Log$

fprintf('    coefs: [%dx%d double]\n', size(f.coefs));
disp(f.classifier)
