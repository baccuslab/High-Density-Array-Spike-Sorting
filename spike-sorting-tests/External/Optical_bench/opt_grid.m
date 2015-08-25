function [opt1] = opt_grid(opt_type,opt_args)
% OPT_GRID - Simple difractive grating
% OPT_TYPE should be 'grid' (hm good and vital argument), OPT_SPEC
% should be a string matrix, see README_OPT for specification.
%   Necessary specification fields: 
% r1, dxdydx, linespmm,
% glass should be name of the glass type, the other 1x3 array
%   Optional specification fields:
% normal, e_slit,
% example OPT_SPEC:
% 'r1          12 0 0'
% 'dxdydz      0 2 1'
% 'linespmm    600'
% 'normal      cos(10*pi/180) sin(10*pi/180) 0'
% 
% Calling:
% [opt1] = opt_grid(opt_type,opt_args)
% 
% See also OPT_APERTURE, OPT_SCREEN, OPT_LENS, OPT_PRISM, OPT_SLIT

% Version: 1.0
% Copyright: Bjorn Gustavsson 20020430


if nargin ~= 2,
  help opt_grid;
  ok = 0;
  return;
end;

opt1 = opt_elem(opt_type);
opt1.glass = 'air'; % 

% defaults
opt1.r(2:3) = 0; % Grid centered on the optical axis
opt1.n = [1 0 0]; % Grid paralell with the optical axis
opt1.e_slits = [0 0 1]; % Grid slits paralell with e_z

% set the necessary ones:
% grid size
ii = opt_findstr(opt_args,'dydz');
opt1.dxdydz = str2num(opt_args(ii,12:end));
% grid position
ii = opt_findstr(opt_args,'r1');
opt1.r = str2num(opt_args(ii,12:end));
% lines per mm
ii = opt_findstr(opt_args,'linespmm');
opt1.lpmm = str2num(opt_args(ii,12:end));

% set the optional ones
ii = opt_findstr(opt_args,'normal');
if length(ii) == 1
  
  opt1.n = str2num(opt_args(ii(1),12:end));
  
end

ii = opt_findstr(opt_args,'e_slit');
if length(ii) == 1
  
  opt1.e_slits = str2num(opt_args(ii(1),12:end));
  
end
