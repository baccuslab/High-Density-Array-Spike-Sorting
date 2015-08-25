function [opt] = opt_prism(opt_type,opt_spec)
% OPT_LENS - Spherical lens single glass.
% OPT_TYPE should be 'prism' (hm good and vital argument), OPT_SPEC
% should be a string matrix, see README_OPT for specification.
%   Necessary specification fields: 
% r1, r2, normal1, normal2, dxdydx, glass
% glass should be name of the glass type, the other 1x3 array
%   Optional specification fields:
% arc, arc1 arc2
% example OPT_SPEC:
% 'normal1     1 0 0'
% 'normal2     cos(30*pi/180) sin(30*pi/180) 0'
% 'r1          12 0 0'
% 'r2          12.5 0 0'
% 'dxdydz1     0 1.5 1.5'
% 'dxdydz2     0 1.7321 1.5'
% 'glass       bk7'
% 'arc         .5'
% 
% Calling:
% [opt] = opt_prism(opt_type,opt_spec)
%
% See also OPT_APERTURE, OPT_GRID, OPT_LENS, OPT_SCREEN, OPT_SLIT

% Version: 1.0
% Copyright: Bjorn Gustavsson 20020430


if nargin ~= 2,
  help opt_prism;
  ok = 0;
  return;
end;

  
opt1 = opt_elem(opt_type);
opt2 = opt_elem(opt_type);

% defaults
opt2.glass = 'air'; % going out into air after the prism
opt1.r(2:3) = 0; % Prism centered on the optical axis
opt1.n = [1 0 0]; % Prism perpendicular to the optical axis
opt2.r(2:3) = 0; % Prism centered on the optical axis
opt2.n = [1 0 0]; % Prism perpendicular to the optical axis
opt1.arc = 0; % no anti reflection coating
opt2.arc = 0; % no anti reflection coating

% set the necessary ones:
% prism size
ii = opt_findstr(opt_spec,'dydz1');
opt1.dxdydz = str2num(opt_spec(ii,12:end));
ii = opt_findstr(opt_spec,'dydz2');
opt2.dxdydz = str2num(opt_spec(ii,12:end));
% prism position
ii = opt_findstr(opt_spec,'r1');
opt1.r = str2num(opt_spec(ii,12:end));
ii = opt_findstr(opt_spec,'r2');
opt2.r = str2num(opt_spec(ii,12:end));
% prism normal1
ii = opt_findstr(opt_spec,'normal1');
opt1.n = str2num(opt_spec(ii,12:end));
% prism normal2
ii = opt_findstr(opt_spec,'normal2');
opt2.n = str2num(opt_spec(ii,12:end));
% prism glass
ii = opt_findstr(opt_spec,'glass');
opt1.glass = strtok(opt_spec(ii,12:end));

% set the optional ones
ii = opt_findstr(opt_spec,'arc');
if length(ii) == 1
  
  opt1.arc = str2num(opt_spec(ii,12:end));
  opt2.arc = str2num(opt_spec(ii,12:end));
  
elseif length(ii)>1
  
  ii = opt_findstr(opt_spec,'arc1');
  opt1.arc = str2num(opt_spec(ii,12:end));
  ii = opt_findstr(opt_spec,'arc2');
  opt2.arc = str2num(opt_spec(ii,12:end));
  
end

opt = [opt1 opt2];
