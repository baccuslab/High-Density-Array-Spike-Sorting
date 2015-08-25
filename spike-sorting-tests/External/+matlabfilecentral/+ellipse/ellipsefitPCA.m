function R = ellipsefitPCA(XY)
%ELLIPSEFIT - form 2D ellipse fit to given x,y data
%in:
%
% XY: Input matrix of 2D coordinates to be fit. XY = [x(:) y(:)]
%
% report: a structure output with the following fields
%    
%    report.PC: principle components
%    report.e: eigenvalues
%    report.xy: center coordinates
%

C = cov(XY);
[V,D] = eig(C);
R.e = diag(D);
R.PC = V;
R.xy = mean(XY);
