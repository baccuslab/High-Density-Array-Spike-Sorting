function rinter = opt_intersection(optelem,rayin)
% OPT_INTERSECTION - Determine the impact point of an optical
% ray RAYIN on an lens system element OPTELEM
%   
% Calling:
% rinter = opt_intersection(optelem,rayin)
% 
% See also OPT_REFRACTION, OPT_TRACE

% Version: 1.0
% Copyright: Bjorn Gustavsson 20020430


if ~isempty(optelem.r)
  l = (optelem.r(1)-rayin.r(1))/rayin.e(1);
  rinter = point_on_line(rayin.r,rayin.e,l);
end

switch optelem.type
 case 'aperture'
  % intersection between a line and a plane
  rinter = rinter+ rayin.e*( ( dot(optelem.r,optelem.n) - dot(rinter,optelem.n) ) / dot(rayin.e,optelem.n) );
  
 case 'grid'
  
  % intersection between a line and a plane
  rinter = rinter+ rayin.e*( ( dot(optelem.r,optelem.n) - dot(rinter,optelem.n) ) / dot(rayin.e,optelem.n) );
  
 case 'lens'
  
  % intersection between a line and a sphere
  [rinter,en] = sphereintersection(rinter,...
				   rayin.e,...
				   optelem.r+optelem.r_o_curv*optelem.n,...
				   optelem.r_o_curv,...
				   optelem.diameter/2,...
				   optelem.n);
 
 case 'prism'
  
  % intersection between a line and a plane
  rinter = rinter+ rayin.e*( ( dot(optelem.r,optelem.n) - dot(rinter,optelem.n) ) / dot(rayin.e,optelem.n) );
  
 case 'screen'
  
  % intersection between a line and a plane
  rinter = rinter+ rayin.e*( ( dot(optelem.r,optelem.n) - dot(rinter,optelem.n) ) / dot(rayin.e,optelem.n) );
  
 case 'slit'
  
  % intersection between a line and a plane
  rinter = rinter+ rayin.e*( ( dot(optelem.r,optelem.n) - dot(rinter,optelem.n) ) / dot(rayin.e,optelem.n) );
  
 otherwise
  
  if 1 < exist('fminu') & exist('fminu') < 7
    l_inter = fminu(optelem.fcn1,0,[],[],rayin.r,rayin.e,'s',optelem.arglist);
  else
    l_inter = fminsearch(optelem.fcn1,0,[],[],rayin.r,rayin.e,'s',optelem.arglist);
  end
  rinter = point_on_line(rayin.r,rayin.e,l_inter);
  
end



function [rinter,en] = sphereintersection(rl,e_in,rsf,curvature,lensradius,en)
% SPHEREINTERSECTION - intersection between a ray and a spherical
% lens surface 
%   

ex = e_in(1);
ey = e_in(2);
ez = e_in(3);
xl = rl(1);
yl = rl(2);
zl = rl(3);
x0 = rsf(1);
y0 = rsf(2);
z0 = rsf(3);
R = curvature;

% I'd like to say that I did this by myself...

l(2) = [ 1/2/(ex^2+ey^2+ez^2)*(2*ez*z0-2*zl*ez-2*yl*ey+2*ex*x0+2*ey*y0-2*xl*ex+2*(2*ex^2*zl*z0+2*ez*z0*ex*x0-ex^2*y0^2-ex^2*z0^2-ex^2*yl^2+ex^2*R^2-ex^2*zl^2-ey^2*xl^2-ey^2*z0^2-ey^2*x0^2+ey^2*R^2-ey^2*zl^2-ez^2*xl^2-ez^2*y0^2-ez^2*x0^2-ez^2*yl^2+ez^2*R^2+2*ex^2*yl*y0+2*ey^2*zl*z0+2*ey^2*xl*x0+2*ez^2*xl*x0+2*ez^2*yl*y0-2*ez*z0*yl*ey+2*ez*z0*ey*y0-2*ez*z0*xl*ex+2*zl*ez*yl*ey-2*zl*ez*ex*x0-2*zl*ez*ey*y0+2*zl*ez*xl*ex-2*yl*ey*ex*x0+2*yl*ey*xl*ex+2*ex*x0*ey*y0-2*ey*y0*xl*ex)^(1/2))];

% ...But who am I trying to fool

l(1) = [ 1/2/(ex^2+ey^2+ez^2)*(2*ez*z0-2*zl*ez-2*yl*ey+2*ex*x0+2*ey*y0-2*xl*ex-2*(2*ex^2*zl*z0+2*ez*z0*ex*x0-ex^2*y0^2-ex^2*z0^2-ex^2*yl^2+ex^2*R^2-ex^2*zl^2-ey^2*xl^2-ey^2*z0^2-ey^2*x0^2+ey^2*R^2-ey^2*zl^2-ez^2*xl^2-ez^2*y0^2-ez^2*x0^2-ez^2*yl^2+ez^2*R^2+2*ex^2*yl*y0+2*ey^2*zl*z0+2*ey^2*xl*x0+2*ez^2*xl*x0+2*ez^2*yl*y0-2*ez*z0*yl*ey+2*ez*z0*ey*y0-2*ez*z0*xl*ex+2*zl*ez*yl*ey-2*zl*ez*ex*x0-2*zl*ez*ey*y0+2*zl*ez*xl*ex-2*yl*ey*ex*x0+2*yl*ey*xl*ex+2*ex*x0*ey*y0-2*ey*y0*xl*ex)^(1/2))];

% Where would I be without the symbolic toolbox?

[qwe,ii] = min(abs(l));
l = l(ii);

rinter = point_on_line(rl,e_in,l);

r1 = rinter - rsf;
r2 = cross(r1,en);

en = [];
if norm(r2)<lensradius
  en = r1/norm(r1);
else
  rinter = [];
end
