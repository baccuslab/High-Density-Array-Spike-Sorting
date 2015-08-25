% Example showing curvature of focus (aberration).

% Set the options for ploting of the optics/raytracing
opt_ops.plotrays=1; % plot the individual rays
opt_ops.plotpaus=0; % do not pause after each  ray intersection
opt_ops.plotRT = 0; % do not ``drawnow'' after each ray intersection

% build the optical system as specified in opt.exmpl
opt_el = opt_build('opt.6.exmpl');


ray = opt_ray;

nr_rays = 3;

for ray_i = 1:nr_rays,
  
  rays(ray_i) = ray;
  rays(ray_i).r = opt_el(1).r+.1*(ray_i-2)*[0 1 0];
  %rays(ray_i).e = [cos(fi_of_i) 0 sin(fi_of_i)];
  rays(ray_i).e = [1 0 0];
  rays(ray_i).color = [0 0 1];% blue - color is line color for
                              % plotting of rays
  
end


r1 = fliplr(1./[-1/6:1/12-eps:1/6]);
r1 = r1(3);

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
  
  r = [-150 3 0];
  dr = [0 -3 0];
  
  for j = 1:5,
    
    opt_el(end).img = 0*opt_el(end).img;
    
    r = r+dr;
    e = -r/norm(r);
    for ray_i = 1:nr_rays,
      
      rays(ray_i).e = e;
      
    end
    
    opt_plotoptics(opt_el,2);
    hold on
    opt_trace(opt_el,rays,opt_ops);
    
  end
  
end
pause(2)
plot([11 12 12 11 11],[-0.2 -0.2 1.1 1.1 -0.2])
disp('push any key')
pause
axis([10.9 12.1 -0.25 1.15])
