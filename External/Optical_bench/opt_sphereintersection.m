function [rinter,en] = opt_sphereintersection(rl,e_in,rsf,curvature,lensradius,en)
% OPT_SPHEREINTERSECTION - intersection between a ray and a spherical lens surface 
%   
% CALLING: 
%  [rinter,en] = opt_sphereintersection(rl,e_in,rsf,curvature,lensradius,en)
% INPUT:
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

