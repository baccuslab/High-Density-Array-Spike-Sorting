% Example showing spherical aberation and coma

% Set the options for ploting of the optics/raytracing
opt_ops.plotrays=0; % plot the individual rays
opt_ops.plotpaus=0; % do not pause after each  ray intersection
opt_ops.plotRT = 0; % do not ``drawnow'' after each ray intersection
% $$$ opt_ops.test = 6; % quick projection calculation

opt_ops2.plotrays=1; % plot the individual rays
opt_ops2.plotpaus=0; % do not pause after each  ray intersection
opt_ops2.plotRT = 0; % do not ``drawnow'' after each ray intersection

ray = opt_ray;

nr_rays = 18;

clrs = jet(9);

for ray_i = 2:nr_rays,
  
  fi_of_i = atan(2/20-ray_i/100);
  
  rays(ray_i) = ray;
  rays(ray_i).r = [-20 2-4*ray_i/20 0];
  %rays(ray_i).e = [cos(fi_of_i) 0 sin(fi_of_i)];
  rays(ray_i).e = [1 0 0];
  rays(ray_i).color = clrs(abs(ray_i-10)+1,:); % color is line color for
                                   % plotting of rays
  
end

% build the optical system as specified in opt.exmpl
%opt_el = opt_build('opt.6.exmpl');
opt_el = opt_build('opt.sca.exmpl');

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
  disp(opt_el(end).r)
  r = [-150 0 0];
  dr = [0 0 0;0 0 16];
  
  for jj = 1:2,
    
    opt_el(end).img = 0*opt_el(end).img;
    
    r = r+dr(jj,:);
    e = -r/norm(r);
    
    %%%subplot(5,4,1+4*(ii-1)+2*(jj-1))
    subplot(5,4,1+4*(ii-1)+1*(jj-1))
    opt_el = opt_project_direction(e,100,ray.wavelength,opt_el,opt_ops);
    opt_plotoptics(opt_el,-2);
    if jj == 1
      axis([-.5 .5 -2 2 -2 2])
      hold on
      opt_trace(opt_el,rays,opt_ops2);
      view(0,90)
      zlnr(1) = round(100*opt_el(2).r_o_curv)/100;
      zlnr(2) = round(100*opt_el(3).r_o_curv)/100;
      ylabel(['r = [',num2str(zlnr(1),'%5g'),', ',num2str(zlnr(2),'%5g'),']'],'fontsize',10)
      if ii == 1
        title('Lens shape','fontsize',16,'FontWeight','bold')
      end
    else
      hold on
      opt_trace(opt_el,rays,opt_ops2);
      view(0,90)
      set(gca,'zticklabel','')
      axis([11 12 -.1 .1 -.1 .1])
      if ii == 1
        title('Rays in focal plane','fontsize',16,'FontWeight','bold')
      end
    end
    
    if ii < 5
      set(gca,'xticklabel','')
    end
    %%%subplot(5,4,2+4*(ii-1)+2*(jj-1))
    subplot(5,4,3+4*(ii-1)+1*(jj-1))
    if jj == 1
      xi = (257-9):(257+9);
      yi = (257-9):(257+9);
    else
      xi = (257-9):(257+9);
      yi = 30:60;
    end
    convK = [.25 .5 .25;.5 1 .5;.25 .5 .25];
    convK = convK/sum(convK);
    
    imagesc(xi,yi,log(conv2(opt_el(end).img(yi,xi),convK,'same')))%,my_cbar
    
    if jj == 1
      sphcax(ii,:) = caxis;
    else
      comacax(ii,:) = caxis;
    end
    if ii == 1
      if jj == 1
        title('Spherical A','fontsize',16,'FontWeight','bold')
      else
        title('Coma','fontsize',16,'FontWeight','bold')
      end
    end
  end
  disp(['ii = ',num2str(ii)])
  drawnow
  
end
