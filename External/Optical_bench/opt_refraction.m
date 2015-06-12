function [en,refraction_loss,n_refr,phase_shift] = opt_refraction(r_int,ray,optelem)
% OPT_REFRACTION - calculation of optical refraction - Snells law.
%   In addition it calculates the fraction of light lost in
%   reflection.
%   
% Calling
% [en,refraction_loss,n_refr,phase_shift] = opt_refraction(r_int,ray,optelem)
% 
% See also OPT_INTERSECTION, OPT_TRACE

% Version: 1.0
% Copyright: Bjorn Gustavsson 20020430

n_refr = 1;
en = ray.e;
refraction_loss = 0;
phase_shift = 0;
fi_out = [];

switch optelem.type
  
 case 'aperture'
  r_int = r_int+ ray.e*( ( dot(optelem.r,optelem.n) - dot(r_int,optelem.n) ) / dot(ray.e,optelem.n) );
  dr = r_int-optelem.r;
  % when falling outside the aperture opening take it away!
  if norm(dr)>optelem.diameter/2
    en = 0*en;
    refraction_loss = 1;
  end
  
 case 'grid'
  
  r_int = r_int+ ( ( dot(optelem.r,optelem.n) - ...
		       dot(r_int,optelem.n) ) / ...
		     dot(ray.e,optelem.n) );
  m = 0;
  n_vec = optelem.n;
  n_vec = n_vec*sign(n_vec(1));
  theta_in = asin(norm(cross(ray.e,n_vec)));
  theta_out = [];
  % Split the input ray to the main maximas.
  while ( abs( sin(theta_in) - m*ray.wavelength*1000*optelem.lpmm ) < 1 )
    
    theta_out = [theta_out,asin( sin(theta_in) - m*ray.wavelength*(1000*optelem.lpmm) )];
    m = m+1;
    
  end
  m = -1;
  % in both directions!
  while ( abs( sin(theta_in) - m*ray.wavelength*1000*optelem.lpmm ) < 1 )
    
    theta_out = [theta_out,asin( sin(theta_in) - m*ray.wavelength*1000*optelem.lpmm)];
    m = m-1;
    
  end
  en = [];
  for i = 1:length(theta_out)
    
    dfi = theta_in-theta_out(i);
    en = [en;opt_rot(ray.e,optelem.e_slits,dfi)];
    refractin_loss = 0;
    
  end
  
 case 'lens'
  % If r_o_curv == inf <=> n_inter == optelem.n
  coc = optelem.r+optelem.n*optelem.r_o_curv;
  r = r_int-coc;
  n_vec = r/norm(r);
  n_vec = n_vec*sign(n_vec(1));
  fi_in = asin(norm(cross(ray.e,n_vec)));
  n_refr = opt_refrindx(optelem.glass,ray.wavelength);
  if ray.n/n_refr*sin(fi_in) < 1
    
    fi_in;
    fi_out = asin(ray.n/n_refr*sin(fi_in));
    dfi = fi_in-fi_out;
    if norm(cross(ray.e,n_vec))
      en = opt_rot(ray.e,cross(ray.e,n_vec),dfi);
    end
    
  end
  
  if ~isempty(fi_out)
    if fi_in+fi_out ~= 0
      % reflexion losses after refraction
      refraction_loss = (sin(dfi)^2/sin(fi_in+fi_out)^2+tan(dfi)^2/tan(fi_in+fi_out)^2)/2;
      
    else
      % and after normal incidence
      refraction_loss = (n_refr - ray.n).^2./(n_refr + ray.n).^2;
      
    end
  else % total relection
    en = 0*en;
    refraction_loss = 1;
  end
  dr = r_int-optelem.r;
  % when falling outside the lens opening take it away!
  if norm(dr)>optelem.diameter/2
    en = 0*en;
    refraction_loss = 1;
  end
  
 case 'flens'
  % Fresnell lens.
  
  coc = optelem.r+optelem.n*optelem.r_o_curv;
  r = r_int-coc;
  n_vec = r/norm(r);
  n_vec = n_vec*sign(n_vec(1));
  fi_in = asin(norm(cross(ray.e,n_vec)));
  n_refr = opt_refrindx(optelem.glass,ray.wavelength);
  if ray.n/n_refr*sin(fi_in) < 1
    
    fi_in;
    fi_out = asin(ray.n/n_refr*sin(fi_in));
    dfi = fi_in-fi_out;
    if norm(cross(ray.e,n_vec))
      en = opt_rot(ray.e,cross(ray.e,n_vec),dfi);
    end
    
  end
  
  if ~isempty(fi_out)
    if fi_in+fi_out ~= 0
      % reflexion losses after refraction
      refraction_loss = (sin(dfi)^2/sin(fi_in+fi_out)^2+tan(dfi)^2/tan(fi_in+fi_out)^2)/2;
      
    else
      % and after normal incidence
      refraction_loss = (n_refr - ray.n).^2./(n_refr + ray.n).^2;
      
    end
  else % total relection
    en = 0*en;
    refraction_loss = 1;
  end
  dr = r_int-optelem.r;
  % when falling outside the lens opening take it away!
  if norm(dr)>optelem.diameter/2
    en = 0*en;
    refraction_loss = 1;
  end
  
 case 'prism'
  
  r_int = r_int+ ray.e*( ( dot(optelem.r,optelem.n) - ...
			   dot(r_int,optelem.n) ) / ...
			 dot(ray.e,optelem.n) );
  n_vec = optelem.n;
  
  n_vec = n_vec*sign(n_vec(1));
  fi_in = asin(norm(cross(ray.e,n_vec)));
  n_refr = opt_refrindx(optelem.glass,ray.wavelength);
  if ray.n/n_refr*sin(fi_in) < 1
    
    fi_out = asin(ray.n/n_refr*sin(fi_in));
    dfi = fi_in-fi_out;
    en = ray.e;
    if norm(cross(ray.e,n_vec))
      en = opt_rot(ray.e,cross(ray.e,n_vec),dfi);
    end
    
  end
  
  if ~isempty(fi_out)
    if fi_in+fi_out ~= 0
      % reflexion losses after refraction
      refraction_loss = (sin(dfi)^2/sin(fi_in+fi_out)^2+tan(dfi)^2/tan(fi_in+fi_out)^2)/2;
      
    else
      % and after normal incidence
      refraction_loss = (n_refr - ray.n).^2./(n_refr + ray.n).^2;
      
    end
  else % total relection
    en = 0*en;
    refraction_loss = 1;
  end
  
 case 'screen'
  
  r_int = r_int+ ray.e*( ( dot(optelem.r,optelem.n) - dot(r_int,optelem.n) ) / dot(ray.e,optelem.n) );
  if ( r_int(2) < optelem.r(2)-optelem.dxdydz(2)/2 | ...
       r_int(2) > optelem.r(2)+optelem.dxdydz(2)/2 | ...
       r_int(3) < optelem.r(3)-optelem.dxdydz(3)/2 | ...
       r_int(3) > optelem.r(3)+optelem.dxdydz(3)/2 )
    en = ray.e;
    refraction_loss = 0;
  else
    en = ray.e;
    refraction_loss = 1;
  end
  
 case 'slit'
  
  r_int = r_int+ ray.e*( ( dot(optelem.r,optelem.n) - dot(r_int,optelem.n) ) / dot(ray.e,optelem.n) );
  if ( r_int(2) < optelem.r(2)-optelem.dxdydz(2)/2 | ...
       r_int(2) > optelem.r(2)+optelem.dxdydz(2)/2 | ...
       r_int(3) < optelem.r(3)-optelem.dxdydz(3)/2 | ...
       r_int(3) > optelem.r(3)+optelem.dxdydz(3)/2 )
    en = ray.e;
    refraction_loss = 0;
  else
    en = 0*ray.e;
    refraction_loss = 1;
  end
  
 otherwise
  
  %fh = inline([optelem.type,'(r_int,s_or_n,arglist)'],'r_int','s_or_n','arglist');
  n_vec = optelem.fcn2(r_int,'n',optelem.arglist);
  fi_in = asin(norm(cross(ray.e,n_vec)));
  n_refr = opt_refrindx(optelem.glass,ray.wavelength);
  if ray.n/n_refr*sin(fi_in) < 1
    
    fi_in;
    fi_out = asin(ray.n/n_refr*sin(fi_in));
    dfi = fi_in-fi_out;
    if norm(cross(ray.e,n_vec))
      en = opt_rot(ray.e,cross(ray.e,n_vec),dfi);
    end
    
  end
  
  if ~isempty(fi_out)
    if fi_in+fi_out ~= 0
      % reflexion losses after refraction
      refraction_loss = (sin(dfi)^2/sin(fi_in+fi_out)^2+tan(dfi)^2/tan(fi_in+fi_out)^2)/2;
      
    else
      % and after normal incidence
      refraction_loss = (n_refr - ray.n).^2./(n_refr + ray.n).^2;
      
    end
  else % total relection
    en = 0*en;
    refraction_loss = 1;
  end
  
end
