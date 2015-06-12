function r2 = opt_r2_of_fnr1_thin(f,n,r1)
% OPT_R2_OF_FNR1_THIN - second radii of curvature given F, N and R1, thin lens
% F - focal width of the lens, N index of refraction of lens
% material, R1 - radius of curvature of the front surface (or
% -1*radius of curvature of the rear surface)  
%
% Calling:
% r2 = opt_r2_of_fnr1_thin(f,n,r1)
%
% Formula from Hecht, Optics

% Copyright B. Gustavsson

r2 = (1./r1-1./(n-1)./f).^-1;
