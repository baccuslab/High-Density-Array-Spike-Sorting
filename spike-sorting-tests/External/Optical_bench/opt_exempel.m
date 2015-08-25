% Example showing the capabilities of Optical_bench:
% Primarily the ray-tracing of an single lens imaging system and
% its point-spread-function

% Set the options for ploting of the optics/raytracing
opt_ops.plotrays=1; % plot the individual rays
opt_ops.plotpaus=0; % do not pause after each  ray intersection
opt_ops.plotRT = 0; % do not ``drawnow'' after each ray intersection

% build the optical system as specified in opt.exmpl
optelements = opt_build('opt.exmpl');

% make a default ray
ray = opt_ray;

%make 100 paralell rays with random displacement in y-z plane
nr_rays = 1000;
for ray_i = 1:nr_rays,
  
  rays(ray_i) = ray;
  rays(ray_i).r = [-3 -1.5+3*rand -1.5+3*rand];
  rays(ray_i).color = rand([3 1]); % color is line color for
                                   % plotting of rays
  
end

clf
subplot(2,2,1)
% plot the optical system
opt_plotoptics(optelements([1:end-4 end]));
hold on
% trace all rays
[opt_screen] = opt_trace(optelements([1:end-4 end]),rays,opt_ops);
view(-15,75)

subplot(2,2,2)
% Last element should be a screen, which stores the image.
imagesc(opt_screen(end).img)
axis([245 265 245 265])

subplot(2,2,3)
ray.r = [-3 -1 0];
opt_plotoptics(optelements);
hold on
[opt_screen] = opt_trace(optelements,ray,opt_ops);
view(-15,75)
subplot(2,2,2)
subplot(2,2,4)
ray.r = [-3 0 0];
[opt_screen] = opt_trace(optelements([end-1 end]),ray,opt_ops);
view(0,90)
