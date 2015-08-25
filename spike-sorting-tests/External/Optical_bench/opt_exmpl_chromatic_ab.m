% Example showing chromatic aberration

% Set the options for ploting of the optics/raytracing
opt_ops.plotrays=1; % plot the individual rays
opt_ops.plotpaus=0; % do not pause after each  ray intersection
opt_ops.plotRT = 0; % do not ``drawnow'' after each ray intersection

% build the optical system as specified in opt.exmpl
optelements = opt_build('opt_single_lens2.exmpl');

% make a default ray
ray = opt_ray;

%make 100 paralell rays with random displacement in y-z plane
nr_rays = 10;

clrs = [0 0 1
        0 1 1
        0 1 0
        .8 .8 0
        1 0 0];
raynr = [1:5 1:5];
wl = [4300:500:6300]*1e-10;

for ray_i = 1:nr_rays,
  
  rays(ray_i) = ray;
  rays(ray_i).r = [-20 1.5-3*floor(ray_i/6) 0];
  rays(ray_i).e = [1 0 0];
  rays(ray_i).color = clrs(raynr(ray_i),:);
  rays(ray_i).wavelength = wl(raynr(ray_i));
  
end

opt_plotoptics(optelements([2 3 end]),2);
hold on
plot3([-2 40],[0 0],[0 0],'k--')
[opt_screen] = opt_trace(optelements,rays,opt_ops);
view(0,90)



%%%to make the figure
% $$$ if 1
% $$$ opt_chromatic_ab_exmpl;axis([-2 25 -1.6 1.6 -1.6 1.6])
% $$$ plot3([21.5 24 24 21.5 21.5],[-.08 -.08 .08 .08 -.08],[1 1 1 1 1],'k','linewidth',1)
% $$$ ax1 = gca;
% $$$ ax2 = axes('position',[.65 .1 .25 .25]);
% $$$ opt_chromatic_ab_exmpl
% $$$ axis([21.51 24.0 -.08 .08 -.079 .08])
% $$$ set(gca,'xtick',[])
% $$$ set(gca,'ytick',[])
% $$$ box on
% $$$ axes(ax1)
% $$$ axis off
% $$$ end
% $$$ print -depsc2 chromatic_abbL.eps
