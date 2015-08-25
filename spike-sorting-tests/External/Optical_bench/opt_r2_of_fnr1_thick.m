function r2 = opt_r2_of_fnr1_thick(f,n,r1,zc)
% r2 = opt_r2_of_fnr1_thick(f,n,r1,zc)
% 
% OPT_R2_OF_FNR1_THICK - second radii of curvature given F, N R1 and ZC, thick lens
% F - focal width of the lens, N index of refraction of lens
% material, R1 - radius of curvature of the front surface (or
% -1*radius of curvature of the rear surface), ZC - centre
% thickness of the lens


fun = inline('(n-1)*(1/r1-1/r2)+zc*(n-1)^2/n/r1/r2-1/f','r2','r1','n','f','zc');

r2 = opt_r2_of_fnr1_thin(f,n,r1);
r2 = fzero(fun,r2,[],r1,n,f,zc);
