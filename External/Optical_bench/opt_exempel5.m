% Example of Spherical abberation

% Set the options for ploting of the optics/raytracing
opt_ops.plotrays=1; % plot the individual rays
opt_ops.plotpaus=0; % do not pause after each  ray intersection
opt_ops.plotRT = 0; % do not ``drawnow'' after each ray intersection

% build the optical system as specified in opt.exmpl
optelements = opt_build('opt_single_lens2.exmpl');

% make a default ray
ray = opt_ray;

%make 100 paralell rays with random displacement in y-z plane
nr_rays = 18;

clrs = jet(9);

r_test2 = [6 12 24 1e6 -12];
r_test3 = -fliplr([6 12 24 1e6 -12]);

  
for ray_i = 2:nr_rays,
  
  fi_of_i = atan(2/20-ray_i/100);
  
  rays(ray_i) = ray;
  rays(ray_i).r = [-20 0 2-4*ray_i/20];
  %rays(ray_i).e = [cos(fi_of_i) 0 sin(fi_of_i)];
  rays(ray_i).e = [1 0 0];
  rays(ray_i).color = clrs(abs(ray_i-10)+1,:); % color is line color for
                                               % plotting of rays
  
end


for o_i = 1:5,
  
  optelements(2).r_o_curv = r_test2(o_i);
  optelements(3).r_o_curv = r_test3(o_i);
  subplot(5,2,2*o_i-1)
  % plot the optical system
  opt_plotoptics(optelements,2);
  hold on
  
  % trace all rays
  [opt_screen] = opt_trace(optelements,rays,opt_ops);
  title(['r1 = ',num2str(optelements(2).r_o_curv),' r2 = ',num2str(optelements(3).r_o_curv)])
  view(-.0,0)
  axis([-1 1 -.2 .2 -2 2])
  if o_i< 5
    set(gca,'xtick',[])
  end
end

for o_i = 1:5,
  
  optelements(2).r_o_curv = r_test2(o_i);
  optelements(3).r_o_curv = r_test3(o_i);
  subplot(5,2,2*o_i)
  opt_plotoptics(optelements);
  hold on
  
  % trace all rays
  [opt_screen] = opt_trace(optelements,rays,opt_ops);
  view(0,0)
  axis([20 24 -.02 .02 -.1 .1])
  
end
