% Example showing astigmatic aberration


% Set the options for ploting of the optics/raytracing
opt_ops.plotrays=1; % plot the individual rays
opt_ops.plotpaus=0; % do not pause after each  ray intersection
opt_ops.plotRT = 0; % do not ``drawnow'' after each ray intersection
opt_ops.plotrayp=0;

% build the optical system as specified in opt.exmpl
opt_el = opt_build('opt.7.exmpl');


ray = opt_ray;

nr_rays = 22;

r = [-150 0 -16*2^.5];
e = -r/norm(r);
  

for ray_i = 1:nr_rays,
  
  rays(ray_i) = ray;
  rays(ray_i).e = e;
  if ray_i <= 11
    rays(ray_i).r = opt_el(1).r+.3*(ray_i-6)*[0 1 0];
    rays(ray_i).color = [0 0 1];% blue - color is line color for
                                % plotting of rays
  else
    rays(ray_i).r = opt_el(1).r+.3*(ray_i-17)*[0 0 1];
    rays(ray_i).color = [1 0 0];% blue - color is line color for
                                % plotting of rays
  end
  rays(ray_i).r = rays(ray_i).r + -3*rays(ray_i).e;
  
end


r1 = fliplr(1./[-1/6:1/12-eps:1/6]);
r1 = r1(3);

ng = opt_refrindx(opt_el(2).glass,ray.wavelength);
f = opt_f_of_r1nr2zc(opt_el(2).r_o_curv,...
                     ng,...
                     opt_el(3).r_o_curv,...
                     norm(opt_el(3).r-opt_el(2).r));


r2 = opt_r2_of_fnr1_thick(f,ng,r1(1),norm(opt_el(3).r-opt_el(2).r));

opt_el(3).r_o_curv = -r1(1);
opt_el(2).r_o_curv = -r2;
F(1) = opt_f_of_r1nr2zc(opt_el(2).r_o_curv,...
                        ng,...
                        opt_el(3).r_o_curv,...
                        norm(opt_el(3).r-opt_el(2).r));

vini = opt_vini(opt_el(2).r_o_curv,opt_el(3).r_o_curv,ng,norm(opt_el(3).r-opt_el(2).r));
opt_el(end).r(1) = opt_el(3).r(1) - vini*opt_el(3).n(1) + F(1)*opt_el(3).n(1);

opt_plotoptics(opt_el([2:3 5:end]),1);
hold on
opt_trace(opt_el([2 3 end]),rays,opt_ops);
alpha = 0.1;
view(-30,20)


keyboard

opt_ops.plotrays=0; % not plot the individual rays

%r_test = [11.582 11.07 10.57 10.82];
r_test = [11.582 10.1 10.2 10.15];
r_test = [11.582 11.2 10.9 11.05];
for jj = 1:length(r_test),
  
  
  opt_el(end).img = 0*opt_el(end).img;
  
  subplot(2,2,jj)
  opt_el(end).r(1) = r_test(jj);
  opt_el = opt_project_direction(e,1,ray.wavelength,opt_el,opt_ops);
  hold off
  cK = [.25 .5 .25;.5 1 .5;.25 .5 .25];
  cK = cK/sum(cK(:));
  imagesc(conv2(opt_el(end).img,cK,'same'))
  drawnow
end
