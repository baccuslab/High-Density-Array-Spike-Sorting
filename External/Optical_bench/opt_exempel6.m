% Example showing Spherical abberation and Coma.


% Set the options for ploting of the optics/raytracing
opt_ops.plotrays=0; % plot the individual rays
opt_ops.plotpaus=0; % do not pause after each  ray intersection
opt_ops.plotRT = 0; % do not ``drawnow'' after each ray intersection

opt_ops2.plotrays=1; % plot the individual rays
opt_ops2.plotpaus=0; % do not pause after each  ray intersection
opt_ops2.plotRT = 0; % do not ``drawnow'' after each ray intersection

ray = opt_ray;

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

% build the optical system as specified in opt.exmpl
opt_el = opt_build('opt.6.exmpl');

r1 = fliplr(1./[-1/12:1/24-eps:1/12]);

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
  
  r = [-150 0 0];
  dr = [0 0 0;0 0 16];
  
  for jj = 1:2,
    
    opt_el(end).img = 0*opt_el(end).img;
    
    r = r+dr(jj,:);
    e = -r/norm(r);
    
    subplot(5,4,1+4*(ii-1)+2*(jj-1))
    opt_el = opt_project_direction(e,100,ray.wavelength,opt_el,opt_ops);
    opt_plotoptics(opt_el);
    if jj == 1
      axis([-1 1 -2 2 -2 2])
      view(0,0)
    else
      hold on
      opt_trace(opt_el,rays,opt_ops2);
    end
    
    title(['r1 = ',num2str(opt_el(2).r_o_curv),' r2 = ',num2str(opt_el(3).r_o_curv)],'fontsize',8)
    set(gca,'xtick',[])
    subplot(5,4,2+4*(ii-1)+2*(jj-1))
    if jj == 1
      xi = (257-9):(257+9);
      yi = (257-9):(257+9);
    else
      xi = (257-9):(257+9);
      yi = 30:60;
    end
    convK = [.25 .5 .25;.5 1 .5;.25 .5 .25];
    convK = convK/sum(convK);
    
    imagesc(xi,yi,log(conv2(opt_el(end).img(yi,xi),convK,'same'))),opt_cbar
    
    if jj == 1
      sphcax(ii,:) = caxis;
    else
      comacax(ii,:) = caxis;
    end
  end
  disp(['ii = ',num2str(ii)])
  
end
