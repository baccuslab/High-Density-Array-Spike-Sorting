% Example of 

figure

% Set the options for ploting of the optics/raytracing
opt_ops.plotrays=1; % plot the individual rays
opt_ops.plotpaus=0; % do not pause after each  ray intersection
opt_ops.plotRT = 0; % do not ``drawnow'' after each ray intersection

% build the optical system as specified in opt.exmpl
optelements = opt_build('opt_single_lens3.exmpl');

% make a default ray
ray = opt_ray;

%make 100 paralell rays with random displacement in y-z plane
nr_rays = 18;

clrs = jet(9);

for ray_i = 2:nr_rays,
  
  fi_of_i = atan(2/20-ray_i/100);
  
  rays(ray_i) = ray;
  rays(ray_i).r = [-20 0 2-4*ray_i/20];
  %rays(ray_i).e = [cos(fi_of_i) 0 sin(fi_of_i)];
  rays(ray_i).e = [1 0 0];
  rays(ray_i).color = clrs(abs(ray_i-10)+1,:); % color is line color for
                                   % plotting of rays
  
end

%subplot(2,2,1)
% plot the optical system
opt_plotoptics(optelements);
hold on
% trace all rays
[opt_screen] = opt_trace(optelements,rays,opt_ops);
view(0,90)

figure 
imagesc(opt_screen(end).img)

% $$$ subplot(2,2,2)
% $$$ % Last element should be a screen, which stores the image.
% $$$ imagesc(opt_screen(end).img)
% $$$ axis([245 265 245 265])
% $$$ 
% $$$ subplot(2,2,3)
% $$$ ray.r = [-3 -1 0];
% $$$ opt_plotoptics(optelements);
% $$$ hold on
% $$$ [opt_screen] = opt_trace(optelements,ray,opt_ops);
% $$$ view(-15,75)
