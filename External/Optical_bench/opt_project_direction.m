function opt_el = opt_project_direction(e,I,wavelengths,opt_el,opt_ops)
% OPT_PROJECT_POINT - Project light from direction E with intensity I throug OPT_EL 
%   
% Calling:
%  opt_el = opt_project_direction(e,I,wavelengths,opt_el,opt_ops)
% 
% Input:
%   E - Direction [3x1] array.
%   I - light intensities [Nx1] array
%   WAVELENGTHS - light wavelength (nm)(m)
%   OPT_EL - optical system [Mx1] array of opt_elem
%   OPT_OPS - struct with options, see OPT_TYPICAL_OPS

if length(opt_el(1).diameter)
  dl = -opt_el(1).diameter/2:opt_el(1).diameter/60:opt_el(1).diameter/2;
else
  dl = -max(opt_el(1).dxdydz)/2:max(opt_el(1).dxdydz)/60:max(opt_el(1).dxdydz)/2;
end

if isfield(opt_ops,'test')
  dl = linspace(min(dl),max(dl),opt_ops.test);
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

for li1 =1:length(dl),
  for li2 = 1:length(dl),
    rinter = rc + e1*dl(li1) + e2*dl(li2);
    x1(li2,li1) = rinter(1);
    y1(li2,li1) = rinter(2);
    z1(li2,li1) = rinter(3);
  end
end


ray = opt_ray;

rays(1:size(x1,1),1:size(x1,2)) = ray;

for ii = 1:size(x1,2),
  
  for j = 1:size(x1,1),
    
    rays(j,ii).r = [x1(j,ii) y1(j,ii) z1(j,ii)];
    rays(j,ii).e = e/norm(e);
    
  end
end

for ij = 1:length(I),
  
  for ii = 1:size(x1,2),
    
    for jj = 1:size(x1,1),
      
      rays(jj,ii).I = I(ij);
      rays(jj,ii).wavelength = wavelengths(ij);
      
    end
  end
  
  [opt_el] = opt_trace(opt_el,rays(:),opt_ops);
  
end
