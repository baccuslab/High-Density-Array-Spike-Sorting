% Example showing pin-hole optics.

theta = [30:2:150]'*pi/180;

r_points(:,2:3) = -[theta.*cos(3*(theta)),theta.*sin(2*(theta)).^2];
r_points(:,1) = -30;
itheta = interp1(cumsum((gradient(r_points(:,2)).^2 + gradient(r_points(:,3)).^2).^.5),...
                 theta,...
                 0:0.075:sum((gradient(r_points(:,2)).^2 + gradient(r_points(:,3)).^2).^.5)); 
clear r_points
r_points(:,2:3) = -[itheta.*cos(3*(itheta));itheta.*sin(2*(itheta)).^2]';
r_points(:,1) = -30;



opt_el = opt_build('opt.pinhole.exmpl');

opt_ops.plotrays=0; % plot the individual rays
opt_ops.plotpaus=0; % do not pause after each  ray intersection
opt_ops.plotRT = 0; % do not ``drawnow'' after each ray intersection


pin_d = [.1 0.035 0.012 0.004];

ray = opt_ray;
k = 2*pi/ray.wavelength;

X = linspace(-opt_el(end).dxdydz(2)/2,opt_el(end).dxdydz(3)/2,512);
Y = linspace(-opt_el(end).dxdydz(2)/2,opt_el(end).dxdydz(3)/2,512);
[X,Y] = meshgrid(X,Y);
I = X;
theta = atan((X.^2+Y.^2).^.5/10);
utheta = unique(theta);                                   

fi = [0:6:360]*pi/180;

for ii = 1:length(pin_d),
  
  %vary the pin hole aperture.
  opt_el(1).diameter = pin_d(ii);
  
  % calculate the diffraction pattern for the apperture
  a = opt_el(1).diameter/200;
  uI = [2*besselj(1,k*a*sin(utheta))./(k*a*sin(utheta))].^2;
  I(:) = interp1(utheta,uI,theta(:));
  % take out the central part - for later convolution.
  if ii < 4
    cK = I(244:269,244:269);
  else
    cK = I(232:281,232:281);
  end
  % additional antialiasing kernel QDF
  cK = conv2(cK,[.25 .5 .25;.5 1 .5;.25 .5 .25],'same');
  
  % show the central diffraction pattern
  subplot(length(pin_d),2,1+2*(ii-1))
  hold off
  imagesc(X(1,:),Y(:,1),I)
  axis([-.1 .1 -.1 .1])
  hold on
  % and the "circle of confusion"
  plot(4/3*opt_el(1).diameter/2*cos(fi),4/3*opt_el(1).diameter/2*sin(fi),'g')
  th(ii) = title(['Aperture diameter: ',num2str(a,'%0.3g'),'(m)'],'fontsize',14);
  if ii < 4
    set(gca,'xticklabel','')
  end
  % zero-set the image
  opt_el(end).img = 0*opt_el(end).img;
  
  subplot(length(pin_d),2,2*(ii))
  for jj = 1:length(r_points),
    
    opt_el = opt_project_point(r_points(jj,:),1,ray.wavelength,opt_el,opt_ops);    
    
  end
  imagesc(X(1,:),Y(:,1),conv2(opt_el(end).img,cK,'same'))
  if ii < 4
    set(gca,'xticklabel','')
  end
  drawnow
  disp(ii)
end
