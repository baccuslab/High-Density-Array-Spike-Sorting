function opt_el = opt_project_point(r,I,wavelengths,opt_el,opt_ops)
% OPT_PROJECT_POINT - Project light from point R with intensity I through optics OPT_EL 
%   
% Calling: 
% opt_el = opt_project_point(r,I,wavelengths,opt_el,opt_ops)
% 
% Input:
%   R - point [3x1] array.
%   I - light intensities [Nx1] array
%   WAVELENGTHS - light wavelength (nm)(m)
%   OPT_EL - optical system [Mx1] array of opt_elem
%   OPT_OPS - struct with options, see OPT_TYPICAL_OPS
% 


if length(opt_el(1).diameter)
  dl = -opt_el(1).diameter/2:opt_el(1).diameter/60:opt_el(1).diameter/2;
else
  dl = -max(opt_el(1).dxdydz)/2:max(opt_el(1).dxdydz)/60:max(opt_el(1).dxdydz)/2;
end
n = opt_el(1).n;
ey = [0 1 0];
rc = opt_el(1).r;

if dot(n,ey)>1-100*eps
  e1 = cross(n,ey);
  e1 = e1/norm(e1);
  e2 = cross(n,e1);
else
  e1 = cross(n,[0 0 1]);
  e1 = e1/norm(e1);
  e2 = cross(n,e1);
end
x1 = nan*zeros([length(dl) length(dl)]);
y1 = x1;
z1 = x1;
%dl(16)
for li1 =1:length(dl),
  for li2 = 1:length(dl),
    rinter = rc + e1*dl(li1) + e2*dl(li2);
    x1(li2,li1) = rinter(1);
    y1(li2,li1) = rinter(2);
    z1(li2,li1) = rinter(3);
  end
end
dA = mean(diff(dl))^2;

% $$$ x1 = x1
% $$$ y1 = y1
% $$$ z1 = z1


ray = opt_ray;

rays(1:size(x1,1),1:size(x1,2)) = ray;

for ii = 1:size(x1,2),
  
  for j = 1:size(x1,1),
    
    rays(j,ii).r = r;
    rays(j,ii).e = [x1(j,ii) y1(j,ii) z1(j,ii)]-r;
    rays(j,ii).e = rays(j,ii).e/norm(rays(j,ii).e);
    
  end % for j = 1:size(x1,
end % for ii = 1:size(x1,

%%%
%%% HAER MAASTE JAG MULTIPLICERA MED cos(theta)*grad(dl)^2
%%% OCH EVENTUELLT MED 1/(4*pi*r^2)
%%% gjort

for ii = 1:length(I),
  
  for iii = 1:size(x1,2),
    
    for jj = 1:size(x1,1),
      
      rays(jj,iii).I = [I(ii) / ...
                      norm([x1(jj,ii) y1(jj,ii) z1(jj,ii)]-r)^2 * ...
                      abs(dot(rays(jj,iii).e,opt_el(1).n)) * ...
                      dA];
      rays(jj,iii).wavelength = wavelengths(ii);
      
    end
  end
  
  [opt_el] = opt_trace(opt_el,rays(:),opt_ops);
  
end
