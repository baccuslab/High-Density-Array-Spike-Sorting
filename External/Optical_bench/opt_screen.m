function [opt1] = opt_screen(opt_type,opt_args)
% OPT_SCREEN - Screen - imaging detector.
% OPT_TYPE should be 'screen' (hm good and vital argument), OPT_SPEC
% should be a string matrix, see README_OPT for specification.
%   Necessary specification fields: 
% r1, dxdydx, imgres
%   Optional specification fields:
% normal
% example OPT_SPEC:
% 'normal      1 0 0'
% 'r1          12 0 0'
% 'dxdydz      0 1.5 1.5'
%   
% Calling:
% [opt1] = opt_screen(opt_type,opt_args)
% 
% See also OPT_LENS, OPT_SCREEN, OPT_GRID, OPT_PRISM, OPT_SLIT

% Version: 1.0
% Copyright: Bjorn Gustavsson 20020430

if nargin ~= 2,
  help opt_screen;
  ok = 0;
  return;
end

opt1 = opt_elem(opt_type);

% defaults
opt1.r(2:3) = 0;   % Screen centered on the optical axis
opt1.n = [1 0 0];   % Screen paralell with the optical axis
opt1.glass = 'air';
% set the necessary ones:
% screen size
ii = opt_findstr(opt_args,'dydz');
opt1.dxdydz = str2num(opt_args(ii,12:end));
% screen position
ii = opt_findstr(opt_args,'r1');
opt1.r = str2num(opt_args(ii,12:end));
ii = opt_findstr(opt_args,'imgres');
opt1.imgres = str2num(opt_args(ii(1),12:end));
opt1.img = zeros(opt1.imgres);

% set the optional ones
ii = opt_findstr(opt_args,'normal');
if length(ii) == 1
  
  opt1.n = str2num(opt_args(ii(1),12:end));
  
end

