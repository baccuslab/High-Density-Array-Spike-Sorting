if 1
clf
grating = linspace(1,-1,40);
grating(3:4:end) = nan;
plot(0*grating,grating,'k','linewidt',2)
grating = linspace(1,-1,40);
hold on
ei = [cos(21*pi/180) sin(21*pi/180)];
eu = [cos(-37*pi/180) sin(-37*pi/180)];
r1_i = point_on_line([0 grating(11)],ei,-1);
r1_u = point_on_line([0 grating(15)],eu,1);
r2_i = point_on_line([0 grating(15)],ei,-1);
r2_u = point_on_line([0 grating(11)],eu,1);
r3_u = point_on_line([0 grating(11)],eu,sin(37*pi/180)*abs(grating(15)-grating(11)));
r3_i = point_on_line([0 grating(11)],ei,-sin(21*pi/180)*abs(grating(15)-grating(11)));
plot([r1_i(1) 0],[r1_i(2) grating(11)],'b','linewidth',2)
plot([r2_i(1) 0],[r2_i(2) grating(15)],'b','linewidth',2)
plot([r1_u(1) 0],[r1_u(2) grating(15)],'b','linewidth',2)
plot([r2_u(1) 0],[r2_u(2) grating(11)],'b','linewidth',2)
plot([r3_u(1) 0],[r3_u(2) grating(11)],'r','linewidth',2)
plot([r3_i(1) 0],[r3_i(2) grating(11)],'r','linewidth',2)
plot([r3_i(1) 0],[r3_i(2) grating(15)],'k--','linewidth',1)
plot([r3_u(1) 0],[r3_u(2) grating(15)],'k--','linewidth',1)
axis equal
axis off
plot([-5*r3_i(1) 0],[1 1]*grating(11),'k--','linewidth',1)
plot([10*r3_i(1) 0],[1 1]*grating(11),'k--','linewidth',1)
gtext('\theta_{in}','fontsize',14)
gtext('\theta_{out}','fontsize',14)
plot([10*r3_i(1) 0],[1 1]*grating(23),'k--','linewidth',1)
plot([10*r3_i(1) 0],[1 1]*grating(27),'k--','linewidth',1)
gtext('d','fontsize',14)
end
