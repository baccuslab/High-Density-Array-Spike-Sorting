function opt_elem = opt_elem(type)
% OPT_ELEM - Default opt_elem structure
% 
% Calling
%  opt_elem = opt_elem(type)
% 
% See also OPT_APERTURE, OPT_GRID, OPT_LENS, OPT_PRISM, OPT_SCREEN, OPT_SLIT

% Version: 1.0
% Copyright: Bjorn Gustavsson 20020430


opt_elem.type = type;

opt_elem.r = []; % location of the optical element
opt_elem.n = []; % normal of the optical element
opt_elem.dxdydz = []; % rectangular size
opt_elem.diameter = []; % Circular size
opt_elem.glass = ''; % Glass type, for apertures and the like -> `air'!
opt_elem.arc = [];   % Anti reflectin coating, currently scalar(*)
opt_elem.imgres = [];% size of image: [sx sy]
opt_elem.img = [];   % image
opt_elem.r_o_curv = []; % spherical curvature of lenses.
opt_elem.lpmm = [];     % lines per milimeter for grids
opt_elem.e_slits = [];  % unit vector || with grid slits
opt_elem.absorption = 0; % absorpption
opt_elem.fcn1 = inline('');
opt_elem.fcn2 = inline('');
opt_elem.arglist = struct;
