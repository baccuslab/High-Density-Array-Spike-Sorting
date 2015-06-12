function ray = opt_ray()
% OPT_RAY - RAY creator for OPT_TOOLS
% Creates a default ray structure  
%   
% CAlling:
% ray = opt_ray()
% 
% See also OPT_REFRACTION, OPT_TRACE

% Version: 1.0
% Copyright: Bjorn Gustavsson 20020430


ray.r = [0 0 0];
ray.e = [1 0 0];
ray.n = 1;
ray.I = 1;
ray.wavelength = 5577e-10;
ray.color = [0 1 0];
ray.absorption = 0;
%ray.glass = 'air'; ???
