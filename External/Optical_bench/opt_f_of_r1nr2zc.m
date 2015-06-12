function f = opt_f_of_r1nr2zc(r1,n,r2,zc)
% OPT_R2_OF_FNR1_THICK - second radii of curvature given F, N R1 and ZC, thick lens
% F - focal width of the lens, N index of refraction of lens
% material, R1 - radius of curvature of the front surface (or
% -1*radius of curvature of the rear surface), ZC - centre
% thickness of the lens
% 
% Calling:
% f = opt_f_of_r1nr2zc(f,n,r1,zc)
% 
% Formula from Hecht, Optics

% Copyright B. Gustavsson

f = 1./((n-1)*(1/r1-1/r2)+zc*(n-1)^2/n/r1/r2);
