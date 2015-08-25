function Hs = opt_plotoptics(opt_el,plotmode)
% OPT_PLOTOPTICS - Plot the optical system.
% 
% Blue - apertures
% Green - lens surfaces
% Red - screens
% Orange - prisms
% Yellow - grids
% 
% Calling:
% h = opt_plotoptics(opt_el,plotmode)
% 
% Input:
%   PLOTMOE - 3 for 3D, 2 for 2D in X-Z plane -2 for 2D in X-Y
%   plane 

% Version: 1.0
% Copyright: Bjorn Gustavsson 20020430

persistent  lens_edge

lens_edge = [];

if nargin == 1
  plotmode = 3;%3-D
end

hstate = ishold;
hold on


for ii = 1:length(opt_el),
  
  switch opt_el(ii).type
   case 'aperture'
    ex = [1 0 0];
    n = opt_el(ii).n;
    dr = opt_el(ii).diameter/2;
    %%% FIXA!!! fixat?
    z = (-1.1*dr):dr/50:(1.1*dr);
    z = sign(z).*((abs(z).^.5));
    z = linspace(0,1,30).^.5;
    z = 1.1*dr*[-z(end:-1:2) z];
    if ii < length(opt_el)
      if length(opt_el(ii+1).diameter)
        z(1) = -abs(opt_el(ii+1).diameter)/2;
        z(end) = abs(opt_el(ii+1).diameter)/2;
      else
        z(1) = -max(opt_el(ii+1).dxdydz)/2;
        z(end) = max(opt_el(ii+1).dxdydz)/2;
      end
    end
    y = z;
    [y,z] = meshgrid(y,z);
    ia = find(y.^2+z.^2<(dr)^2);
    x = zeros(size(y));
    r = [x(:),y(:),z(:)];
    if asin(norm(cross(ex,n)))~=0
      for jj = 1:length(r)
	r(jj,:) = opt_rot(r(jj,:),cross(ex,n),asin(norm(cross(ex,n))));
      end
    end
    x(:) = r(:,1)+opt_el(ii).r(1);
    y(:) = r(:,2)+opt_el(ii).r(2);
    z(:) = r(:,3)+opt_el(ii).r(3);
    y(ia) = nan;
    
    if plotmode == 3
      Hs(ii).h = surf(x,y,z,zeros(size(z)));
    else
      ex = [1 0 0];
      n = opt_el(ii).n;
      dr = opt_el(ii).diameter/2;
      %%% FIXA!!! fixat?
      z = [-1.1*dr -dr 0  dr 1.1*dr];
      [y,z] = meshgrid(z,z);
      x = zeros(size(y));
      r = [x(:),y(:),z(:)];
      if asin(norm(cross(ex,n)))~=0
        for jj = 1:length(r)
          r(jj,:) = opt_rot(r(jj,:),cross(ex,n),asin(norm(cross(ex,n))));
        end
      end
      x(:) = r(:,1)+opt_el(ii).r(1);
      y(:) = r(:,2)+opt_el(ii).r(2);
      z(:) = r(:,3)+opt_el(ii).r(3);
      y(3,3) = nan;
% $$$       disp([x(3,:);y(3,:);z(3,:)])
% $$$       y(:,3) = nan;
      Hs(ii).h = plot3(x(3,:),y(3,:),z(3,:),x(:,3),y(:,3),z(:,3),'b');
% $$$       keyboard
    end
    shading flat
    
   case 'grid'
    ex = [1 0 0];
    n = opt_el(ii).n;
    dz = opt_el(ii).dxdydz(3);
    dy = opt_el(ii).dxdydz(2);

    z = (-dz/2):(dz/2):(dz/2);
    y = (-dy/2):(dy/2):(dy/2);
    [y,z] = meshgrid(y,z);
    x = zeros(size(y));
    r = [x(:),y(:),z(:)];
    if asin(norm(cross(ex,n)))~=0
      for jj = 1:length(r)
	r(jj,:) = opt_rot(r(jj,:),cross(ex,n),asin(norm(cross(ex,n))));
      end
    end
    x(:) = r(:,1)+opt_el(ii).r(1);
    y(:) = r(:,2)+opt_el(ii).r(2);
    z(:) = r(:,3)+opt_el(ii).r(3);
    if plotmode == 3
      Hs(ii).h = surf(x,y,z,.6*ones(size(z)));shading flat
      plot3([opt_el(ii).r(1)  opt_el(ii).r(1)+opt_el(ii).e_slits(1)],...
            [opt_el(ii).r(2)  opt_el(ii).r(2)+opt_el(ii).e_slits(2)],...
            [opt_el(ii).r(3)  opt_el(ii).r(3)+opt_el(ii).e_slits(3)],'k','linewidth',6)
    else
      Hs(ii).h = plot3([opt_el(ii).r(1)  opt_el(ii).r(1)+opt_el(ii).e_slits(1)],...
                      [opt_el(ii).r(2)  opt_el(ii).r(2)+opt_el(ii).e_slits(2)],...
                      [opt_el(ii).r(3)  opt_el(ii).r(3)+opt_el(ii).e_slits(3)],'k','linewidth',6);
    end
    
   case 'lens'
    
    ex = [1 0 0];
    n = opt_el(ii).n;
    r = opt_el(ii).diameter/2;
    
    R = 0:r/20:r;
    fi = [0:10:360]*pi/180;
    [fi,R] = meshgrid(fi,R);
    
    y = R.*cos(fi);
    z = R.*sin(fi);
    
% $$$     x = -(sign(opt_el(ii).r_o_curv)*(opt_el(ii).r_o_curv^2- ...
% $$$                                     y.^2-z.^2).^.5 + ...
% $$$           -sign(opt_el(ii).r_o_curv)*(opt_el(ii).r_o_curv^2-4*r^2).^.5);
    x = -sign(opt_el(ii).r_o_curv)*(opt_el(ii).r_o_curv^2- ...
                                    y.^2-z.^2).^.5;
          %-sign(opt_el(ii).r_o_curv)*(opt_el(ii).r_o_curv^2-4*r^2).^.5);
% $$$     if sign(opt_el(ii).r_o_curv) > 0
% $$$       x = x-min(x(:));
% $$$     else
% $$$       x = x-max(x(:));
% $$$     end
    x = x-(x(1,19));
    
    r = [x(:),y(:),z(:)];
    if asin(norm(cross(ex,n)))~=0
      [ex,n,cross(ex,n),asin(norm(cross(ex,n)))]
      for jj = ib,
	r(jj,:) = opt_rot(r(jj,:),cross(ex,n),asin(norm(cross(ex,n))));
      end
    end
    x(:) = r(:,1)+opt_el(ii).r(1);
    y(:) = r(:,2)+opt_el(ii).r(2);
    z(:) = r(:,3)+opt_el(ii).r(3);
    x(11,19)
    %surf(x,y,z,0.5*ones(size(z))),shading flat
% $$$     surf(x,y,z,10*(x-min(x(:)))),shading flat
% $$$     plot3(x,y,z),shading flat
    if ~strcmp(opt_el(ii).glass,'air')
      lens_edge = [x([end],[1 10 19 28]),y([end],[1 10 19 28]),z([end],[1 10 19 28])];
      lens_edge = reshape(lens_edge,[4 3]);
    else
      lens_edge2 = [x([end],[1 10 19 28]),y([end],[1 10 19 28]),z([end],[1 10 19 28])];
      lens_edge2 = reshape(lens_edge2,[4 3]);
      switch plotmode
       case 3
        plot3([lens_edge(:,1),lens_edge2(:,1)]',[lens_edge(:,2),lens_edge2(:,2)]',[lens_edge(:,3),lens_edge2(:,3)]','b')
       case 2
        plot3([lens_edge([2 4],1),lens_edge2([2 4],1)]',[lens_edge([2 4],2),lens_edge2([2 4],2)]',[lens_edge([2 4],3),lens_edge2([2 4],3)]','b')
       case -2
        plot3([lens_edge([1 3],1),lens_edge2([1 3],1)]',[lens_edge([1 3],2),lens_edge2([1 3],2)]',[lens_edge([1 3],3),lens_edge2([1 3],3)]','b')
       otherwise
        plot3([lens_edge([1 3],1),lens_edge2([1 3],1)]',[lens_edge([1 3],2),lens_edge2([1 3],2)]',[lens_edge([1 3],3),lens_edge2([1 3],3)]','b')
      end
    end
    
    half_end = floor(size(x)/2);
    if plotmode == 3
      Hs(ii).h = surf(x,y,z,0.5*ones(size(z)));shading flat
      plot3(x(1:half_end(1),[1 19]),y(1:half_end(1),[1 19]),z(1:half_end(1),[1 19]),'b');
      plot3(x(1:half_end(1),[10 28]),y(1:half_end(1),[10 28]),z(1:half_end(1),[10 28]),'b');
      plot3(x(end,:),y(end,:),z(end,:),'b');
    elseif plotmode == 2;
      h = plot3(x(1:end,[1 19]),y(1:end,[1 19]),z(1:end,[1 19]),'b');
      hh = plot3(x(1:end,[10 28]),y(1:end,[10 28]),z(1:end,[10 28]),'b');
      hhh = plot3(x(end,:),y(end,:),z(end,:),'b');
      Hs(ii).h = [h;hh;hhh];
    elseif plotmode == -2;
      %x(11,[10 19 28])
      h = plot3(x(1:end,[10 28]),z(1:end,[10 28]),y(1:end,[10 28]),'b');
      Hs(ii).h = h;
      %keyboard
    else
      h = plot3(x(1:half_end(1),[1 19]),y(1:half_end(1),[1 19]),z(1:half_end(1),[1 19]),'b');
      hh = plot3(x(1:half_end(1),[10 28]),y(1:half_end(1),[10 28]),z(1:half_end(1),[10 28]),'b');
      hhh = plot3(x(end,:),y(end,:),z(end,:),'b');
      Hs(ii).h = [h;hh;hhh];
    end
    
   case 'lens_cartesian'
    
    ex = [1 0 0];
    n = opt_el(ii).n;
    dz = opt_el(ii).diameter;
    
    z = (-(dz/2)):dz/30:(dz/2);
    y = z;
    [y,z] = meshgrid(y,z);
    ia = find(y.^2+z.^2>=(dz/2)^2);
    y(ia) = nan;
    ib = find(y.^2+z.^2<(dz/2)^2);
    x = zeros(size(y));
    
    x(ib) = -(sign(opt_el(ii).r_o_curv)*(opt_el(ii).r_o_curv^2- ...
				      y(ib).^2-z(ib).^2).^.5 + ...
	    -sign(opt_el(ii).r_o_curv)*(opt_el(ii).r_o_curv^2-dz^2).^.5);
    if sign(opt_el(ii).r_o_curv) > 0
      x(ib) = x(ib)-min(x(ib));
    else
      x(ib) = x(ib)-max(x(ib));
    end
    r = [x(:),y(:),z(:)];
    if asin(norm(cross(ex,n)))~=0
      [ex,n,cross(ex,n),asin(norm(cross(ex,n)))]
      for jj = ib,
	r(jj,:) = opt_rot(r(jj,:),cross(ex,n),asin(norm(cross(ex,n))));
      end
    end
    x(:) = r(:,1)+opt_el(ii).r(1);
    y(:) = r(:,2)+opt_el(ii).r(2);
    z(:) = r(:,3)+opt_el(ii).r(3);
    if plotmode == 3
      Hs(ii).h = surf(x,y,z,0.5*ones(size(z)));shading flat
    else
      Hs(ii).h = plot3(x(:,16),y(:,16),z(:,16),'b');
    end
    
   case 'lens_newer'
    
    dl = -opt_el(ii).diameter/2:opt_el(ii).diameter/30:opt_el(ii).diameter/2;
    n = opt_el(ii).n;
    ey = [0 1 0];
    rc = opt_el(ii).r;
    
    if dot(n,ey)>1-100*eps
      e1 = cross(n,ey);
      e1 = e1/norm(e1);
      e2 = cross(n,e1);
    else
      e1 = cross(n,[0 0 1]);
      e1 = e1/norm(e1);
      e2 = cross(n,e1);
    end
    x = nan*zeros([length(dl) length(dl)]);
    y = x;
    z = x;
    %dl(16)
    for li1 =1:length(dl),
      for li2 = 1:length(dl),
        rinter = rc + e1*dl(li1) + e2*dl(li2);
        rinter = opt_sphereintersection(rinter,...
                                        opt_el(ii).n,...
                                        opt_el(ii).r+opt_el(ii).r_o_curv*opt_el(ii).n,...
                                        opt_el(ii).r_o_curv,...
                                        opt_el(ii).diameter/2,...
                                        opt_el(ii).n);
        if length(rinter)
          x(li2,li1) = rinter(1);
          y(li2,li1) = rinter(2);
          z(li2,li1) = rinter(3);
        end
      end
    end
    %surf(x,y,z,0.5*ones(size(z))),shading flat
    if plotmode == 3
      Hs(ii).h = surf(x,y,z,0.5*ones(size(z)));shading flat
    else
      Hs(ii).h = plot3(x(16,:),y(16,:),z(16,:),x(:,16),y(:,16),z(:,16),'b');
    end
    
    
   case 'prism'
    ex = [1 0 0];
    n = opt_el(ii).n;
    dz = opt_el(ii).dxdydz(3);
    dy = opt_el(ii).dxdydz(2);
    z = (-dz/2):(dz/2):(dz/2);
    y = (-dy/2):(dy/2):(dy/2);
    [y,z] = meshgrid(y,z);
    x = zeros(size(y));
    r = [x(:),y(:),z(:)];
    if asin(norm(cross(ex,n)))~=0
      for jj = 1:length(r)
	r(jj,:) = opt_rot(r(jj,:),cross(ex,n),asin(norm(cross(ex,n))));
      end
    end
    x(:) = r(:,1)+opt_el(ii).r(1);
    y(:) = r(:,2)+opt_el(ii).r(2);
    z(:) = r(:,3)+opt_el(ii).r(3);
    if plotmode == 3
      Hs(ii).h = surf(x,y,z,0.7*ones(size(z)));shading flat
    else
      Hs(ii).h = plot3(x,y,z,'r');shading flat
    end
   case 'screen'
    ex = [1 0 0];
    n = opt_el(ii).n;
    dz = opt_el(ii).dxdydz(3);
    dy = opt_el(ii).dxdydz(2);
    z = (-dz/2):(dz/2):(dz/2);
    y = (-dy/2):(dy/2):(dy/2);
    [y,z] = meshgrid(y,z);
    x = zeros(size(y));
    r = [x(:),y(:),z(:)];
    if asin(norm(cross(ex,n)))~=0
      for jj = 1:length(r)
	r(jj,:) = opt_rot(r(jj,:),cross(ex,n),asin(norm(cross(ex,n))));
      end
    end
    x(:) = r(:,1)+opt_el(ii).r(1);
    y(:) = r(:,2)+opt_el(ii).r(2);
    z(:) = r(:,3)+opt_el(ii).r(3);
    if plotmode == 3
      Hs(ii).h = surf(x,y,z,ones(size(z)));shading flat
      plot3(x(1,:),...
            y(1,:),...
            z(1,:),'b',...
            x(:,1),...
            y(:,1),...
            z(:,1),'b',...
            x(end,:),...
            y(end,:),...
            z(end,:),'b',...
            x(:,end),...
            y(:,end),...
            z(:,end),'b','linewidth',1);
    else
      Hs(ii).h = plot3(x(1,:),...
                      y(1,:),...
                      z(1,:),'b',...
                      x(:,1),...
                      y(:,1),...
                      z(:,1),'b',...
                      x(end,:),...
                      y(end,:),...
                      z(end,:),'b',...
                      x(:,end),...
                      y(:,end),...
                      z(:,end),'b','linewidth',1);
    end
    
   case 'slit'
    ex = [1 0 0];
    n = opt_el(ii).n;
    dz = opt_el(ii).dxdydz(3)+2;
    dy = opt_el(ii).dxdydz(2)+2;
    z = (-dz/2):(dz/2):(dz/2);
    y = (-dy/2):(dy/2):(dy/2);
    [y,z] = meshgrid(y,z);
    x = zeros(size(y));
    r = [x(:),y(:),z(:)];
    if asin(norm(cross(ex,n)))~=0
      for jj = 1:length(r)
	r(jj,:) = opt_rot(r(jj,:),cross(ex,n),asin(norm(cross(ex,n))));
      end
    end
    x(:) = r(:,1)+opt_el(ii).r(1);
    y(:) = r(:,2)+opt_el(ii).r(2);
    z(:) = r(:,3)+opt_el(ii).r(3);
    size(x)
    if plotmode == 3
      Hs(ii).h = surf(x,y,z,.9*ones(size(z)));shading flat
    else
      Hs(ii).h = plot3([x(1,1) x(1,end) x(end,end) x(end,1) x(1,1)],...
                      [y(1,1) y(1,end) y(end,end) y(end,1) y(1,1)],...
                      [z(1,1) z(1,end) z(end,end) z(end,1) z(1,1)]);
    end
    ex = [1 0 0];
    n = opt_el(ii).n;
    dz = opt_el(ii).dxdydz(3);
    dy = opt_el(ii).dxdydz(2);
    z = (-dz/2):(dz/2):(dz/2);
    y = (-dy/2):(dy/2):(dy/2);
    [y,z] = meshgrid(y,z);
    x = zeros(size(y));
    r = [x(:),y(:),z(:)];
    if asin(norm(cross(ex,n)))~=0
      for jj = 1:length(r)
	r(jj,:) = opt_rot(r(jj,:),cross(ex,n),asin(norm(cross(ex,n))));
      end
    end
    x(:) = r(:,1)+opt_el(ii).r(1);
    y(:) = r(:,2)+opt_el(ii).r(2);
    z(:) = r(:,3)+opt_el(ii).r(3);
    plot3(x,y,z,'k'),shading flat
    plot3(x',y',z','k'),shading flat
    
   otherwise
    feval(opt_el(ii).fcn2,[1 1 1],'p',opt_el(ii).arglist);
  end

end

if ~hstate
  
  hold off
  
end
