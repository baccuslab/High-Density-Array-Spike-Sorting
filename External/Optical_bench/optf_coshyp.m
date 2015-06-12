function out_arg = optf_coshyp(r,s_or_n_or_p,arglist)
% OPTF_HYPERBOLC defines a cosh lens surface
% Called with optf_hyperbolic(R,'s',arglist) the function should
% return 0 (zero) when R is on the lens surface and monotonically
% growing scalars when R deviates. When R is outside the lens area
% but on the analytical extension of the lens surface the function
% should return a scalar smaller than -2eps or larger than 2eps.
% Called with  optf_hyperbolic(R,'n',arglist) the functnion should
% return the normal of the lens surface. If you want the surface to
% appear in a plot of the optics the function should respond with a
% proper surface plot when called with
% optf_hyperbolic(R,'p',arglist) ARGLIST will be a struct as
% produced by OPT_FCN.
%
% Calling:
% function out_arg = optf_hyperbolic(r,s_or_n,arglist)


% parsing of the arglist structure.
r0 = arglist.r0;                  % The lenssurface reference point
lens_rot = arglist.rotmat;        % Rotation matrix of the lens
f = arglist.focals;               % The shape parameters
lens_diameter = arglist.diameter; % Diameter of the lens. 

r = lens_rot'*(r-r0)';

switch s_or_n_or_p
 case 's'
  % calculate the some distance between r and the surface
  % equation for the surface f(x,y,z) = 0
  out_arg = abs(f(3)*(cosh((r(2).^2/f(1).^2+r(3).^2/f(2).^2).^.5)-1)-r(1));
  if r(1)^2+r(2)^2 > lens_diameter^2/4;
    out_arg = out_arg+1;
  end
 case 'n'
  % calculate the surface normal at point r
  % The surface normal is -grad(f)
  out_arg = [-1 ...
             f(3)*r(2)/f(1)^2./((r(2).^2/f(1).^2+r(3).^2/f(2).^2).^.5)*sinh((r(2).^2/f(1).^2+r(3).^2/f(2).^2).^.5) ...
             f(3)*r(3)/f(2)^2./((r(2).^2/f(1).^2+r(3).^2/f(2).^2).^.5)*sinh((r(2).^2/f(1).^2+r(3).^2/f(2).^2).^.5)];
  out_arg = -out_arg/norm(out_arg);
 case 'p'
  % plot the surface
  y = -lens_diameter/2:lens_diameter/40:lens_diameter/2;
  z = -lens_diameter/2:lens_diameter/40:lens_diameter/2;
  [y,z] = meshgrid(y,z);
  x = f(3)*(cosh((y.^2/f(1).^2+z.^2/f(2)^2).^.5)-1);
  R = [x(:) y(:) z(:)];
  ii = find(R(:,2).^2+R(:,3).^2>lens_diameter^2/4);
  R(ii,1) = nan;
  R(ii,2) = nan;
  R(ii,3) = nan;
  R = lens_rot'*R';
  
  x(:) = R(1,:)+r0(1);
  y(:) = R(2,:)+r0(2);
  z(:) = R(3,:)+r0(3);
  %outarg needs to be set, and what more apropriate then the handle?
  out_arg = surf(x,y,z,ones(size(z))*rand(1));
  shading faceted;
 otherwise
  error('Value of s_or_n_or_p (',s_or_n_or_p,') out of range (''s'',''n'',''p'')')
end
