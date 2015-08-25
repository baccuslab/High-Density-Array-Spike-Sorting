function vini = opt_vini(r1,r2,ng,zc)
% OPT_VINI - Rear nodal point distance from rear vertex of thick lens
%   
% Calling:
%  vini = opt_vini(r1,r2,ng,zc)
% 
% Formula from Hecht Optics

% Copyright B Gustavsson 20050404

vini = zc*r2/(ng*(r2-r1)+zc*(ng-1));
