function [optelements] = opt_trace(optelements,rays,opt_ops)
% OPT_TRACE - ray tracing through optical systems.
% This function handles aberationand and attenuation in an optical
% lens system. Currently suported optical elements are Lenses,
% apertures, grids prisms (without internal reflection), slits, and
% transmission grids (only main maximas). Currently unsupported:
% Mirrors, beam splitters, diffraction.
% 
% Calling:
% [optelements] = opt_trace(optelements,rays,opt_ops)
% 

% Version: 1.09
% Copyright: Bjorn Gustavsson 20020430

% trace all rays
for ii = 1:length(rays)
  
  if length(rays(ii).r) & length(rays(ii).e) & rays(ii).I>0
    
    % determine where on the next optical element the ray intersects
    ri_int = opt_intersection(optelements(1),rays(ii));
    if ~isempty(ri_int)
      % calculate:
      % the line-of-sight after refraction
      % the reflected intensity (unpolarised)
      % and the refractive index after the surface
      [en,refrloss,n_refr] = opt_refraction(ri_int,rays(ii),optelements(1));
      % calculate the intensity after absorption and reflection
      I = rays(ii).I*(1-refrloss)*exp(-rays(ii).absorption* ...
                                      (norm(rays(ii).r-ri_int)));
      if isfield(rays(ii),'phase')
        % if the ray caries phase information - update
        rays.phase = ( rays.phase + ...
                       2*pi / r_tmp(ii).n/r_tmp(ii).wavelength * norm(r_tmp(ii).r-ri_int) );
      end
      for jj = 1:size(en,1),
        
        r_tmp(jj) = rays(ii);
        r_tmp(jj).r = ri_int;
        r_tmp(jj).e = en(jj,:);
        r_tmp(jj).I = I/size(en,2);
        r_tmp(jj).absorption = opt_absorption(optelements(1).glass,r_tmp(jj).wavelength);
        r_tmp(jj).n = n_refr;
        
      end
      
      % if ray fall onto screen - add intensity to image.
      if strcmp('screen',optelements(1).type)
        % where in the image plane
        imgindx1 = round(1+( ri_int(2) - (optelements(1).r(2)-optelements(1).dxdydz(2)/2) )/optelements(1).dxdydz(2)*(optelements(1).imgres(1)));
        imgindx2 = round(1+( ri_int(3) - (optelements(1).r(3)-optelements(1).dxdydz(3)/2) )/optelements(1).dxdydz(3)*(optelements(1).imgres(2)));
        % if ray fall onto screen - add intensity to image.
        if ( imgindx1 > 0 & ...
             imgindx2 > 0 & ...
             imgindx1 <= optelements(1).imgres(1) & ...
             imgindx2 <= optelements(1).imgres(2) )
          if isfield(rays(ii),'phase')
            % seems as infinite coherence length to me...
            optelements(1).img((imgindx2),(imgindx1)) = ( optelements(1).img(imgindx2,imgindx1) +...
                                                          rays(ii).I*exp(i*rays(ii).phase) );
          else
            optelements(1).img((imgindx2),(imgindx1)) = ( optelements(1).img(imgindx2,imgindx1) +...
                                                          rays(ii).I );
          end
        end
      end
      % plot the raytrace...
      if isfield(opt_ops,'plotrays') & opt_ops.plotrays==1 & prod(size(en))
        
        plot3([rays(ii).r(1) r_tmp(1).r(1)],...
              [rays(ii).r(2) r_tmp(1).r(2)],...
              [rays(ii).r(3) r_tmp(1).r(3)],'-','color',rays(ii).color)
        if strcmp('screen',optelements(1).type)
% $$$         plot3([rays(ii).r(1) r_tmp(1).r(1)],...
% $$$               [rays(ii).r(2) r_tmp(1).r(2)],...
% $$$               [rays(ii).r(3) r_tmp(1).r(3)],'.','color',rays(ii).color)
% $$$         plot3([rays(ii).r(1)],...
% $$$               [rays(ii).r(2)],...
% $$$               [rays(ii).r(3)],'.','color',rays(ii).color)
          plot3([r_tmp(1).r(1)],...
                [r_tmp(1).r(2)],...
                [r_tmp(1).r(3)],'.','color',rays(ii).color)
        end % if strcmp('screen',
        grid on
        hold on
        % ...and pause...
        if isfield(opt_ops,'plotpaus') & opt_ops.plotpaus
          disp('push any button')
          pause
        end
        % or at least do it real time.
        if isfield(opt_ops,'plotRT') & opt_ops.plotRT
          drawnow
        end
      end
      if isfield(opt_ops,'plotrayp') & opt_ops.plotrayp==1 & prod(size(en))
        
        if strcmp('screen',optelements(1).type)
          plot3([rays(ii).r(1) r_tmp(1).r(1)],...
                [rays(ii).r(2) r_tmp(1).r(2)],...
                [rays(ii).r(3) r_tmp(1).r(3)],'.','color',rays(ii).color)
        end % if strcmp('screen',
      end % if isfield(opt_ops,
      
      % Recursion - sorry about that.
      if length(optelements) > 1 & sum([r_tmp(:).I]) > 0
        [opt_es] = opt_trace(optelements(2:end),r_tmp,opt_ops);
        optelements = [optelements(1) opt_es];
      end
    end
    
  end % if length(rays(ii).
  
end
