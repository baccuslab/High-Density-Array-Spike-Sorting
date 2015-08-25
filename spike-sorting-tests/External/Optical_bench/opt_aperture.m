function [opt1] = opt_aperture(opt_type,opt_args)
% OPT_APERTURE Circular aperture, iris
% OPT_TYPE should be 'aperture' (hm good and vital argument), OPT_SPEC
% should be a string matrix, see README_OPT for specification.
%   Necessary specification fields: 
% r1, r2, normal1, normal2, diameter,
%   Optional specification fields:
% normal
% example OPT_SPEC:
% 'normal     1 0 0'
% 'r1         12 0 0'
% 'diameter   1.5'
% 
% Calling:
% [opt1] = opt_aperture(opt_type,opt_args)
% 
% See also OPT_LENS, OPT_SCREEN, OPT_GRID, OPT_PRISM, OPT_SLIT

% Version: 1.0
% Copyright: Bjorn Gustavsson 20020430

if nargin ~= 2,
  help opt_aperture;
  ok = 0;
  return;
end;


opt1 = opt_elem(opt_type);

% defaults
opt1.r(2:3) = 0;   % Aperture centered on the optical axis
opt1.n = [1 0 0];   % Aperture perpendicular with the optical axis
opt1.arc = 0; % no anti reflection coating
opt1.glass = 'air'; % aperture ``glass''

% set the necessary ones:
% aperture diameter
ii = opt_findstr(opt_args,'diameter');
opt1.diameter = str2num(opt_args(ii,12:end));
% apperture position
ii = opt_findstr(opt_args,'r1');
opt1.r = str2num(opt_args(ii,12:end));

% set the optional ones
ii = opt_findstr(opt_args,'normal');
if length(ii) == 1
  
  opt1.n = str2num(opt_args(ii(1),12:end));
  
end
