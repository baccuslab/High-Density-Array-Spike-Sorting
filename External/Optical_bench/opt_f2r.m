function [r] = opt_f2r(f,n)
% OPT_F2R Thin lens focal width to lens curvature radius.
% F - focal width, N refractive index.
% 
% Calling:
% [opt] = opt_f2r(opt_type,opt_spec)
% 

if nargin ~= 2,
  help opt_f2r;
  ok = 0;
  return;
end;

% Scool-book formula for thin lens:
% 1/F = (n_g-1)(1/r1-1/r2)
% We take it as a symetrical lens (r2 = -r1):
%1/F = (n_g-1)(2/r1)
r = f*(n-1)*2;
