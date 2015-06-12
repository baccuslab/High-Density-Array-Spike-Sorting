function [opt1] = opt_slit(opt_type,opt_args)
% OPT_SLIT - optical slit.
% OPT_TYPE should be 'slit' (hm good and vital argument), OPT_SPEC
% should be a string matrix, see README_OPT for specification.
%   Necessary specification fields: 
% r1, normal, dxdydx,
% glass should be name of the glass type, the other 1x3 array
%   Optional specification fields:
% e_slit
% example OPT_SPEC:
% 'normal     1 0 0'
% 'r1          12 0 0'
% 'dxdydz      0 1.5 60e-6'
% 
% Calling:
% [opt1] = opt_slit(opt_type,opt_args)
% 
% See also OPT_APERTURE, OPT_LENS, OPT_PRISM, OPT_GRID, OPT_SCREEN

% Version: 1.0
% Copyright: Bjorn Gustavsson 20020430


if nargin ~= 2,
  help opt_slit;
  ok = 0;
  return;
end;


opt1 = opt_elem(opt_type);

% defaults
opt1.r(2:3) = 0;   % Slit centered on the optical axis
opt1.n = [1 0 0];   % Slit perpendicular to the optical axis

% set the necessary ones:
% slit size
ii = opt_findstr(opt_args,'dydz');
opt1.dxdydz = str2num(opt_args(ii,12:end));
% slit position
ii = opt_findstr(opt_args,'r1');
opt1.r = str2num(opt_args(ii,12:end));

% set the optional ones
ii = opt_findstr(opt_args,'normal');
if length(ii) == 1
  
  opt1.n = str2num(opt_args(ii(1),12:end));
  
end

% set the optional ones
ii = opt_findstr(opt_args,'e_slit');
if length(ii) == 1
  
  opt1.e_slits = str2num(opt_args(ii(1),12:end));
  
end
