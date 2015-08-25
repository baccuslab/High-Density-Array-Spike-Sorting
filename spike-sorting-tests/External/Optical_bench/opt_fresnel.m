function [opt] = opt_fresnel(opt_type,opt_spec)
% OPT_LENS - Spherical lens single glass.
% OPT_TYPE should be 'lens' (hm good and vital argument), OPT_SPEC
% should be a string matrix, see README_OPT for specification.
%   Necessary specification fields: 
% r1, r2, curvature1, curvature2, dxdydx, glass, diameter
% glass should be name of the glass type, the other 1x3 array
%   Optional specification fields:
% arc, arc1 arc2, normal
% example OPT_SPEC:
% 'curvature1  12'
% 'curvature2  -15'
% 'r1          12 0 0'
% 'r2          12.5 0 0'
% 'diameter    4'
% 'glass       bk7'
% 'arc         .5'
%   
% Calling
% [opt] = opt_lens(opt_type,opt_spec)
% 
% See also OPT_APERTURE, OPT_SCREEN, OPT_GRID, OPT_PRISM, OPT_SLIT

% Version: 1.0
% Copyright: Bjorn Gustavsson 20020430

if nargin ~= 2,
  help opt_lens;
  ok = 0;
  return;
end;

opt1 = opt_elem('lens');
opt2 = opt_elem('lens');

% defaults
opt2.glass = 'air'; % going out into air after the lens
opt1.r(2:3) = 0;   % Lens centered on the optical axis
opt2.r(2:3) = 0;   % Lens centered on the optical axis
opt1.n = [1 0 0];   % Lens perpendicular to the optical axis
opt2.n = [1 0 0];   % Lens perpendicular to the optical axis
opt1.arc = 0; % no anti reflection coating
opt2.arc = 0; % no anti reflection coating

% set the necessary ones:
% lens diameter
ii = opt_findstr(opt_spec,'diameter');
opt1.diameter = str2num(opt_spec(ii,12:end));
opt2.diameter = str2num(opt_spec(ii,12:end));
% lens glass 
ii = opt_findstr(opt_spec,'glass');
opt1.glass = strtok(opt_spec(ii,12:end));
% lens front radius of curvature
ii = opt_findstr(opt_spec,'curvature1');
opt1.r_o_curv = str2num(opt_spec(ii,12:end));
% lens back radius of curvature
ii = opt_findstr(opt_spec,'curvature2');
opt2.r_o_curv = str2num(opt_spec(ii,12:end));
% lens front position
ii = opt_findstr(opt_spec,'r1');
opt1.r = str2num(opt_spec(ii,12:end));
% lens back position
ii = opt_findstr(opt_spec,'r2');
opt2.r = str2num(opt_spec(ii,12:end));

% set the optional ones
ii = opt_findstr(opt_spec,'normal');
if length(ii) == 1
  
  opt1.n = str2num(opt_spec(ii,12:end));
  opt2.n = str2num(opt_spec(ii,12:end));
  
elseif length(ii)>1
  
  ii = opt_findstr(opt_spec,'normal1');
  opt1.n = str2num(opt_spec(ii,12:end));
  ii = opt_findstr(opt_spec,'normal2');
  opt2.n = str2num(opt_spec(ii,12:end));
  
end

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

ii = opt_findstr(opt_args,'linespmm');
opt1.lpmm = str2num(opt_args(ii,12:end));

opt = [opt1 opt2];
