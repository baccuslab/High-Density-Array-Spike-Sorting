% Example showing coma (aberration)

% Set the options for ploting of the optics/raytracing
opt_ops.plotrays=1; % plot the individual rays
opt_ops.plotrayp=0; % plot the individual rays intersection points
                    % as points
opt_ops.plotpaus=0; % do not pause after each  ray intersection
opt_ops.plotRT = 0; % do not ``drawnow'' after each ray intersection

opt_ops2.plotrays=0; % plot the individual rays
opt_ops2.plotrayp=1; % plot the individual rays intersection points
                    % as points
opt_ops2.plotpaus=0; % do not pause after each  ray intersection
opt_ops2.plotRT = 0; % do not ``drawnow'' after each ray intersection

ray = opt_ray;

% build the optical system as specified in opt.6.exmpl
opt_el = opt_build('opt.6.exmpl');

r = [-150 0 -16];

nr_rays = 3;

clrs = jet(7);
clrs = [1 0 0;0 1 0;0 0 1];
fi_ray = [0:360/16:359]*pi/180;
r_frac = [.2 .5 .8];
for ray_i = 1:nr_rays,
  
  for ray_j = 1:length(fi_ray)
    
    
    rays(ray_i,ray_j) = ray;
    yr = cos(fi_ray(ray_j)) * opt_el(1).diameter/2*r_frac(4-ray_i);
    zr = sin(fi_ray(ray_j)) * opt_el(1).diameter/2*r_frac(4-ray_i);
    rays(ray_i,ray_j).r = [opt_el(1).r(1) yr zr];
    %rays(ray_i,ray_j).e = [1 0 0];%opt_el(1).r - r;
    rays(ray_i,ray_j).e = opt_el(1).r - r;
    rays(ray_i,ray_j).e = rays(ray_i,ray_j).e/norm(rays(ray_i,ray_j).e);
    %rays(ray_i,ray_j).r = rays(ray_i,ray_j).r - 4*rays(ray_i,ray_j).e;
    
    rays(ray_i,ray_j).color = clrs(ray_i,:); % color is line color for
                                   % plotting of rays
  end
end


r1 = fliplr(1./[-1/6:1/12-eps:1/6]);
r1 = r1([1 3 4]);
ng = opt_refrindx(opt_el(2).glass,ray.wavelength);
f = opt_f_of_r1nr2zc(opt_el(2).r_o_curv,...
                     ng,...
                     opt_el(3).r_o_curv,...
                     norm(opt_el(3).r-opt_el(2).r));

for ii = 1:length(r1),
  
  r2 = opt_r2_of_fnr1_thick(f,ng,r1(ii),norm(opt_el(3).r-opt_el(2).r));
  
  opt_el(3).r_o_curv = -r1(ii);
  opt_el(2).r_o_curv = -r2;
  F(ii) = opt_f_of_r1nr2zc(opt_el(2).r_o_curv,...
                          ng,...
                          opt_el(3).r_o_curv,...
                          norm(opt_el(3).r-opt_el(2).r));
  
  vini = opt_vini(opt_el(2).r_o_curv,opt_el(3).r_o_curv,ng,norm(opt_el(3).r-opt_el(2).r));
  opt_el(end).r = opt_el(3).r - vini*opt_el(3).n + F(ii)*opt_el(3).n;
  
  
  opt_el(end).img = 0*opt_el(end).img;
  
  subplot(3,5,1+5*(ii-1))
  %opt_el = opt_project_direction(e,100,ray.wavelength,opt_el,opt_ops);
  hs1 = opt_plotoptics(opt_el(2:3),2);
  Hs1(ii,1:length(hs1)) = hs1;
  axis([-.3 .3 -1.5 1.5 -1.5 1.5])
  view(0,0)
  %axis off
  %set(gca,'xticklabel','')
  subplot(3,3,2+3*(ii-1))
  hold on
  hs2 = opt_plotoptics(opt_el([2 3 end]),1);
  Hs2(ii,1:length(hs2)) = hs2;
  
  opt_screen = opt_trace(opt_el([2 3 end]),rays(end:-1:1),opt_ops);
  title(['r1 = ',num2str(opt_el(2).r_o_curv),' r2 = ',num2str(opt_el(3).r_o_curv)],'fontsize',12)
  set(gca,'xticklabel','','yticklabel','','zticklabel','')
  view(90-35,15)
  alpha(.2)
  axis([-3 12 -2 2 -1.7 1.7])
% $$$   plot3(opt_el(end).r(1)+1000*eps,0,1.25,'b.','markersize',14)
  pax = get(gca,'position');
  set(gca,'position',pax+[-.12 0 .2 0])
  
  subplot(3,4,4+4*(ii-1))
  hs2 = opt_plotoptics(opt_el,2);
  Hs3(ii,1:length(hs2)) = hs2;
  hold on
  
  opt_screen = opt_trace(opt_el,rays(:),opt_ops2);
  axis image
  axis([11 12 -.2 .2 1.05 1.45])
  view(90,0)
  %alpha(.2)
% $$$   xi = (257-9):(257+9);
% $$$   yi = 30:60;
% $$$   imagesc(xi,yi,opt_screen(end).img(yi,xi))
   disp(['ii = ',num2str(ii)])
  
end

% $$$ set(Hs1(1,2).h(2),'color',[1 1 1],'linewidth',0)
% $$$ set(Hs1(2,1).h(2),'color',[1 1 1],'linewidth',0)
% $$$ set(Hs1(3,1).h(2),'color',[1 1 1],'linewidth',0)
% $$$ set(Hs1(3,2).h(2),'color',[1 1 1],'linewidth',0)
set(Hs1(1,2).h(2),'color',[1 1 1])
set(Hs1(2,1).h(2),'color',[1 1 1])
set(Hs1(3,1).h(2),'color',[1 1 1])
set(Hs1(3,2).h(2),'color',[1 1 1])
set(Hs1(1,2).h(1),'color',[1 1 1])
set(Hs1(2,1).h(1),'color',[1 1 1])
set(Hs1(3,1).h(1),'color',[1 1 1])
set(Hs1(3,2).h(1),'color',[1 1 1])
