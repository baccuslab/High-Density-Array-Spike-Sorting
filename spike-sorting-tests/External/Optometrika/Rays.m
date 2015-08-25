classdef Rays
    % RAYS Implements a ray bundle
    % Note that for easy copying Rays doesn't inherit from handle
    %
    % Member functions:
    %
    % r = Rays( n, geometry, r, dir, D, pattern, glass, wavelength, color ) - object constructor
    % INPUT:
    %   n - number of rays in the bundle
    %   geometry - For geometry 'collimated', r defines rays origins, while dir - 
    % their direction. For geometry 'source', r defines position 
    % of the point source, and dir - direction along which rays propagate 
    %   r - 1x3 bundle source position vector
    %   dir - 1x3 bundle direction vector
    %   D - diameter of the ray bundle (at distance 1 if geometry = 'source' )
    %   pattern - (optional) pattern of rays within the bundle: 'linear', 'hexagonal'
    % or 'random', hexagonal by default
    %   glass - (optional) material through which rays propagate, 'air' by
    % default
    %   wavelength - (optional) wavelength of the ray bundle, meters, 557.7
    % nm by default
    %   color - (optional) 1 x 3 vector defining the color with which to draw
    % the ray, [ 0 1 0 ] (green) by default
    %
    % OUTPUT:
    %   r - ray bundle object
    %
    % r.draw( scale ) - draws the ray bundle r in the current axes as arrows
    % INPUT:
    %   scale - (optional) the arrow length, 1 by default
    %
    % [ rays_out, nrms ] = r.intersection( surf ) - finds intersection
    % of the ray bundle r with a surface
    % INPUT:
    %   surf - the surface
    % OUTPUT:
    %   rays_out - rays_out.r has the intersection points, the remaining
    % structure parameters might be wrong
    %   nrms - 1x3 normal vectors to the surface at the intersection points
    %
    % rays_out = r.interaction( surf ) - finishes forming the outcoming
    % ray bundle, calculates correct directions and intensities
    % INPUT:
    %   surf - surface
    % OUTPUT:
    %   rays_out - outcoming ray bundle
    %
    % r = r.append( r1 ) - appends to bundle r bundle r1
    % INPUT:
    %   bundle - the appended ray bundle
    % OUTPUT:
    %   r - the resulting ray bundle
    % 
    % r = r.truncate() - truncate all rays with zero intensity from the bundle r.
    %
    % [ av, dv, nrays ] = r.stat() - return statistics on the bundle r
    % OUTPUT:
    %   av - 1x3 vector of the mean bundle position
    %   dv - standard deviation of the ray positions in the bundle
    %   nrays - number of rays with non-zero intensity in the bundle
    %
    % [ x0, cv, ax, ang, nrays ] = r.stat_ellipse() - fit a circumscribing ellipse
    % to the bundle r in the YZ plane
    % OUTPUT:
    %   x0 - 1x3 vector of the ellipse center
    %   cv - bundle covariance matrix
    %   ax - 1x2 vector of the ellipse half-axes lengths
    %   ang - angle of rotation of the ellipse from the longer axis being
    %   oriented along the Y axis.
    %   nrays - number of rays with non-zero intensity in the bundle
    %
    % r2 = dist2rays( p ) - returns squared distances from point p to all rays
    % INPUT:
    %   p - 1x3 vector
    % OUTPUT:
    %   r2 - nrays x 1 vector of squared distances
    %
    % f = r.focal_point() - find a focal point of the bundle. The focal
    % point is defined as the location in space where the bundle cross-section is the tightest.
    % OUTPUT:
    %    f - 1x3 vector for the focal point 
    %
    % Copyright: Yury Petrov, Oculus VR, 01/2014
    %
    
    properties
        r = [] % a matrix of ray starting positions
        n = [] % a matrix of ray directions
        w = 5577e-10  % a vector of ray wavelengths
        I = 1         % a vector of ray intensities
        nrefr = []    % a vector of current refractive indices
        att = 0       % a vector of ray attenuations
        color = [ 0 1 0 ];   % color to draw the bundle rays
        cnt = 1;      % number of rays in the bundle
    end
    
    methods
        function self = Rays( cnt, geometry, pos, dir, diameter, rflag, glass, wavelength, acolor ) % constructor of ray bundles
            % Constructs a ray bundle comprising 'cnt' rays. For geometry 
            % 'collimated', 'pos' defines rays origins, while 'dir' - 
            % their direction. For geometry 'source', 'pos' defines position 
            % of the point source, 'dir' - direction along which rays form a 
            % linear, hexagonal, or random pattern (specified by 'rflag') of 
            % the size specified by 'diameter' at distance 1.
            
            if nargin == 0 % used to allocate arrays of Rays
                return;
            end
                
            if nargin > 8
                self.color = acolor;
            end
            if nargin > 7
                self.w = wavelength;
            end
            if nargin < 7
                glass = 'air';
            end
            if nargin < 6
                rflag = 'hexagonal'; % hexagonal lattice of rays
            end
            if nargin < 5
                diameter = 1;
            end
            if nargin < 4
                dir = [ 1 0 0 ];
            end
            if nargin < 3
                pos = [ 0 0 0 ];
            end
            if nargin < 2
                geometry = 'collimated';
            end
                    
            if strcmp( rflag, 'linear' ) % extend along y-axis
                p( :, 1 ) = linspace( -diameter/2, diameter/2, cnt ); % all rays starting from the center
                p( :, 2 ) = 0;
            elseif strcmp( rflag, 'random' )
                cnt1 = round( cnt * 4 / pi );
                p( :, 1 ) = diameter * ( rand( cnt1, 1 ) - 0.5 ); % horizontal positions
                p( :, 2 ) = diameter * ( rand( cnt1, 1 ) - 0.5 ); % vertical positions
                p( p( :, 1 ).^2 + p( :, 2 ).^2 > diameter^2 / 4, : ) = []; % leave rays only within the diameter
            elseif strcmp( rflag, 'hexagonal' )
                % find the closest hexagonal number to cnt
                cnt1 = round( cnt * 2 * sqrt(3) / pi );
                tmp = (-3 + sqrt( 9 - 12 * ( 1 - cnt1 ) ) ) / 6;
                cn( 1 ) = floor( tmp );
                cn( 2 ) = ceil(  tmp );
                totn = 1 + 3 * cn .* ( 1 + cn );
                [ ~, i ] = min( abs( totn - cnt1 ) );
                cn = cn( i );
                % generate hexagonal grid
                p = [];
                for i = cn : -1 : -cn % loop over rows starting from the top
                    nr = 2 * cn + 1 - abs( i ); % number in a row
                    hn = floor( nr / 2 );
                    if rem( nr, 2 ) == 1
                        x = ( -hn : hn )';
                    else
                        x = ( -hn : hn - 1 )' + 1/2;
                    end
                    p = [ p; [ x, i * sqrt( 3 ) / 2 * ones( nr, 1 ) ] ]; % add new pin locations
                end
                if cn > 0
                    p = p * diameter / 2 / cn * 2 / sqrt( 3 ); % circubscribe the hexagon by an inward circle
                end
                if cn > 2 % cut away corners of the hexagon
                    p( p( :, 1 ).^2 + p( :, 2 ).^2 > ( diameter / ( 4 * cn ) )^2 * 4 / 3 + diameter^2 / 4, : ) = [];
                end
            elseif strcmp( rflag, 'pentile' )
                % generate a pentile grid
                dim = round( sqrt( cnt / 8 ) ); % number of cells in each dimension
                p = zeros( dim^2 * 8, 2 );
                p( 1, : ) = [ 0   0 ]; % green origin
                p( 2, : ) = [ .5   0 ]; % green right-bottom
                p( 3, : ) = [ 0   .5 ]; % green left-top
                p( 4, : ) = [ .5  .5 ]; % green right-top
                p( 5, : ) = [ .25 .25 ]; % red left-bottom
                p( 6, : ) = [ .75 .75 ]; % red right-top
                p( 7, : ) = [ .25 .75 ]; % blue left-top
                p( 8, : ) = [ .75 .25 ]; % blue right-bottom
                self.w( 1:4, 1 ) = 5300e-10;  % green wavelength of the OLED display
                self.color( 1:4, : ) = repmat( [ 0 1 0 ], 4, 1 );
                self.w( 5:6, 1 ) = 6200e-10;  % red wavelength of the OLED display
                self.color( 5:6, : ) = repmat( [ 1 0 0 ], 2, 1 );
                self.w( 7:8, 1 ) = 4580e-10;  % blue wavelength of the OLED display
                self.color( 7:8, : ) = repmat( [ 0 0 1 ], 2, 1 );
                
                for i = 0 : dim - 1
                    for j = 1 : dim
                        ind = 8 * ( i * dim + j );
                        p( ind + 1 : ind + 8, 1 ) = p( 1:8, 1 ) + j - ceil( dim/2 );
                        p( ind + 1 : ind + 8, 2 ) = p( 1:8, 2 ) + i - floor( dim/2 );
                        self.w( ind + 1 : ind + 8, 1 ) = self.w( 1:8 );
                        self.color( ind + 1 : ind + 8, : ) = self.color( 1:8, : );
                    end
                end
                p = p( 9 : end, : ); % remove the original cell
                self.w = self.w( 9 : end );
                self.color = self.color( 9 : end, : );
                p = p * diameter / 2 / max( p( :, 1 ) ); % scale the ray positions
                ind = p( :, 1 ).^2 + p( :, 2 ).^2 > diameter^2 / 4;
                p( ind, : ) = []; % leave rays only within the diameter
                self.w( ind ) = [];
                self.color( ind, : ) = [];
            else
                error( [ 'Ray arrangement flag ' rflag ' is not defined!' ] );
            end
            self.cnt = size( p, 1 );
            
            p = [ zeros( self.cnt, 1 ) p ]; % add x-positions
            pos = repmat( pos, self.cnt, 1 );
            % normalize direction
            dir = dir ./ norm( dir );
            if strcmp( geometry, 'collimated' ) % parallel rays
                % distribute over the area
                self.r = pos + p;
                dir = repmat( dir, self.cnt, 1 );
                self.n = dir;
            elseif strcmp( geometry, 'source' ) || strcmp( geometry, 'source-Lambert' ) % assume p array at dir, source at pos.
                ex = [ 1 0 0 ];
                ax = cross( ex, dir );
                if norm( ax ) ~= 0
                    p = rodrigues_rot( p, ax, asin( norm( ax ) ) );
                end
                self.r = pos;
                self.n = p + repmat( dir, self.cnt, 1 );
            else
                error( [ 'Source geometry' source ' is not defined!' ] );
            end
            % normalize directions
            self.n = self.n ./ repmat( sqrt( sum( self.n.^2, 2 ) ), 1, 3 );
           
            if ~strcmp( rflag, 'pentile' )
                self.w = repmat( self.w, self.cnt, 1 );
                self.color = repmat( self.color, self.cnt, 1 );
            end
            self.nrefr = refrindx( self.w, glass );
            self.I = ones( self.cnt, 1 );
            if strcmp( geometry, 'source-Lambert' )
                self.I = self.I .* self.n( :, 1 ); % Lambertian source: I proportional to cos wrt source surface normal assumed to be [ 1 0 0 ]
            end
            
            self.att = ones( self.cnt, 1 );
        end
            
        function draw( self, scale )
            if nargin == 0 || isempty( scale )
                scale = 1;
            end
            vis = self.I ~= 0;
            [ unique_colors, ~, ic ] = unique( self.color, 'rows' );
            nrms = scale * self.n;
            for i = 1 : size( unique_colors, 1 )
                cvis = vis & ( ic == i );
                quiver3( self.r( cvis, 1 ), self.r( cvis, 2 ), self.r( cvis, 3 ), ...
                         nrms( cvis, 1 ),   nrms( cvis, 2 ),   nrms( cvis, 3 ), ...
                        0, 'Color', unique_colors( i, : ), 'ShowArrowHead', 'off' );
            end
        end
         
        function [ rays_out, nrms ] = intersection( self, surf )
            % instantiate Rays object
            rays_out = self; % copy incoming rays

            switch class( surf )

                case { 'Aperture', 'Plane', 'Screen' } % intersection with a plane                        
                    % distance to the plane along the ray
                    d = dot( repmat( surf.n, self.cnt, 1 ), repmat( surf.r, self.cnt, 1 ) - self.r, 2 ) ./ ...
                        dot( self.n, repmat( surf.n, self.cnt, 1 ), 2 );
                    
                    % calculate intersection vectors and normals
                    rinter = self.r + repmat( d, 1, 3 ) .* self.n;
                    nrms = repmat( surf.n, self.cnt, 1 );
                    
                    % bring surface to the default position
                    rtr = rinter - repmat( surf.r, self.cnt, 1 );
                    if surf.rotang ~= 0
                        rtr = rodrigues_rot( rtr, surf.rotax, -surf.rotang ); % rotate rays to the default plane orientation
                    end
              
                    if isa( surf, 'Aperture' )
                        % mark intersections outside the aperture diameter as zero intensity
                        ind = sum( rtr( :, 2:3 ).^2, 2 ) > ( surf.D / 2 )^2;
                        rays_out.I( ind ) = 0;
                        % leave the starting points unchanGed for the rays passing through the aperture
                        rays_out.r( ind, : ) = rinter( ind, : );
                        % rays_out.r = rinter;
                    elseif isa( surf, 'Plane' )
                        rays_out.r = rinter;
                        % mark intersections outside the plane boudnaries as zero intensity
                        if ~( surf.w == 0 && surf.h == 0 )
                            ind =  rtr( :, 2 ) < -surf.w/2 | rtr( :, 2 ) > surf.w/2 | ...
                                   rtr( :, 3 ) < -surf.h/2 | rtr( :, 3 ) > surf.h/2;
                        elseif surf.D ~= 0
                            ind = sum( rtr( :, 2:3 ).^2, 2 ) > ( surf.D / 2 )^2;
                        end
                        rays_out.I( ind ) = 0;
                    elseif isa( surf, 'Screen' ) % calculate retinal image
                        rays_out.r = rinter;
                        % mark intersections outside the plane boudnaries as zero intensity
                        ind =  rtr( :, 2 ) < -surf.w/2 | rtr( :, 2 ) > surf.w/2 | ...
                               rtr( :, 3 ) < -surf.h/2 | rtr( :, 3 ) > surf.h/2; 
                        rays_out.I( ind ) = 0;
                        rays_out.r( ind, : ) = Inf; % this is not to draw rays that missed the screen
                        % form image
                         surf.image = hist2( rtr( :, 2 ), rtr( :, 3 ), self.I, ...
                                             linspace( -surf.w/2, surf.w/2, surf.wbins ), ...
                                             linspace( -surf.h/2, surf.h/2, surf.hbins ) );
                         surf.image = flipud( surf.image ); % to get from matrix to image form
                         surf.image = fliplr( surf.image ); % because y-axis points to the left
                    end                
                    
                case { 'GeneralLens', 'FresnelLens', 'Lens', 'Retina' } % intersection with a conical surface of rotation
                    % intersection between rays and the surface, also returns surface normals at the intersections
                    
                    % transform rays into the lens surface RF
                    r_in = self.r - repmat( surf.r, self.cnt, 1 ); % shift to RF with surface origin at [ 0 0 ]
                    if surf.rotang ~= 0 % rotate so that the surface axis is along [1 0 0]
                        r_in = rodrigues_rot( r_in, surf.rotax, -surf.rotang ); % rotate rays to the default surface orientation
                        e = rodrigues_rot( self.n, surf.rotax, -surf.rotang );
                    else
                        e = self.n;                        
                    end
                   
                    if isa( surf , 'GeneralLens' )
                        % minimize a measure of distance between a ray point and the surface
                        rinter = ones( self.cnt, 3 ); % init intersection vectors
                        if exist( 'fminunc', 'file' ) % requires optimization toolbox
                            options = optimoptions( 'fminunc', 'Algorithm', 'quasi-newton', 'Display', 'off', 'Diagnostics', 'off' );
                            parfor i = 1 : self.cnt % run parallel computing
                                rinter( i, : ) = r_in( i, : ) + e( i, : ) * fminunc( @dist2, 0, options, r_in( i, : ), e( i, : ), surf );
                            end
                        else % no optimization toolbox
                            options = optimoptions( 'fminunc', 'MaxFunEvals', 2000, 'Display', 'off', 'Diagnostics', 'off' );
                            parfor i = 1 : self.cnt  % run parallel computing
                                rinter( i, : ) = r_in( i, : ) + e( i, : ) * fminsearch( @dist2, 0, options, r_in( i, : ), e( i, : ), surf );
                            end
                        end
                        % get surface normals at the intersection points
                        en = surf.funch( rinter( :, 2 ), rinter( :, 3 ), surf.funca, 1 ); 
                    
                    elseif isa( surf , 'FresnelLens' ) % Fresnel lens
                        % find rays intersections with each Fresnel cone
                        rays_out.I = 0 * rays_out.I;
                        rinter = Inf * ones( self.cnt, 3 ); % init intersection vectors
                        rings = ones( self.cnt, 1 );        % ring indices of the rays
                        for i = 1 : surf.ncones
                            the = surf.the( i );
                            cs2 = cos( the )^2;
                            c2 = e( :, 1 ).^2 - cs2; % d^2 coefficient
                            % find the cone's vertex coordinates and cone slice limits
                            minsag = surf.sag( i );
                            if the == pi/2 % frontoparallel ring
                                dm = r_in( :, 1 ) - surf.sag( i );
                                dp = Inf;
                            else
                                if i == 1
                                    vx = surf.sag( i );
                                    maxsag = surf.sag( i ) + surf.rad( i ) / tan( the );
                                else
                                    vx = surf.sag( i ) - surf.rad( i - 1 ) / tan( the );
                                    maxsag = surf.sag( i ) + ( surf.rad( i ) - surf.rad( i - 1 ) ) / tan( the );
                                end
                                if minsag == maxsag
                                    maxsag = minsag + realmin; % make the two tiny different to avoid rays disappearing at frontoparallel surfaces
                                elseif minsag > maxsag
                                    tmp = minsag;
                                    minsag = maxsag;
                                    maxsag = tmp;
                                end
                                v = [ vx 0 0 ]; % cone vertex coordinates
                                dt = r_in - repmat( v, size( r_in, 1 ), 1 ); % rays in the reference frame of the cone vertex at [ 0 0 0 ]
                                c1 = 2 * ( e( :, 1 ) .* dt( :, 1 ) - cs2 * dot( e, dt, 2 ) ); % d^1 coefficient
                                c0 = dt( :, 1 ).^2 - cs2 * sum( dt.^2, 2 ); % d^0 coefficient
                                
                                % solve c2 * d^2 + c1 * d + c0 = 0 for d
                                D =  c1.^2 - 4 * c2 .* c0;
                                D( D < -1e-12 ) = Inf; % suppress negative values (no intersection with the half-cone) to avoid imaginary values
                                D( D < 0 ) = 0; % remove round off negative Ds
                                D = sqrt( D );
                                dm = 0.5 * ( -c1 - D ) ./ c2;
                                dp = 0.5 * ( -c1 + D ) ./ c2;
                            end
                            
                            % form the intersection vector
                            rm = r_in + repmat( dm, 1, 3 ) .* e;
                            rp = r_in + repmat( dp, 1, 3 ) .* e;
                            % tmp = [ rm; rp ]; figure, scatter3( tmp(:,1), tmp(:,2), tmp(:,3) ), axis equal vis3d
                            
                            % save intersections inside the cone ring for either of the two d values
                            r2 = sum( rm( :, 2 : 3 ).^2, 2 );
                            if i == 1
                                in = r2 <= surf.rad( i )^2 & rm( :, 1 ) > minsag & rm( :, 1 ) <= maxsag;
                            else
                                in = r2 <= surf.rad( i )^2 & r2 > surf.rad( i - 1 )^2 & ...
                                     rm( :, 1 ) > minsag & rm( :, 1 ) <= maxsag;
                            end
                            rays_out.I( in, : ) = 1;
                            rinter( in, : ) = rm( in, : );
                            rings( in, 1 ) = i; % memorize ring indices for the intersecting rays
                            
                            r2 = sum( rp( :, 2 : 3 ).^2, 2 );
                            if i == 1
                                in = r2 <= surf.rad( i )^2 & rp( :, 1 ) > minsag & rp( :, 1 ) <= maxsag;
                            else
                                in = r2 <= surf.rad( i )^2 & r2 > surf.rad( i - 1 )^2 & ...
                                     rp( :, 1 ) > minsag & rp( :, 1 ) <= maxsag;
                            end
                            rays_out.I( in, : ) = 1;
                            rinter( in, : ) = rp( in, : );                            
                            rings( in, 1 ) = i; % memorize ring indices for the intersecting rays                            
                       end                      
                       % find normals
                       c = cos( surf.the( rings ) );
                       s = sqrt( 1 - c.^2 );
                       th = atan2( rinter( :, 3 ), rinter( :, 2 ) ); % rotation angle to bring r into XZ plane
                       en = [ s, -c .* cos( th ), -c .* sin( th ) ]; % make normal sign positive wrt ray
                       % Correct for rays hitting at the very center of the
                       % lens, where the refraction is underined. Assume
                       % that the refraction does not happen
                       central_rays = sum( rinter( :, 2:3 ).^2, 2 ) < 1e-20; % rays hitting the very center of the Fresnel lens
                       if sum( central_rays ) > 0
                           en( central_rays, : ) = self.n( central_rays, : );
                       end
                       % figure, scatter3( rinter( :, 1 ), rinter( :, 2 ), rinter( :, 3 ) ), hold on, quiver3( rinter( :, 1 ), rinter( :, 2 ), rinter( :, 3 ), en( :, 1 ), en( :, 2 ), en( :, 3 ), 5 ), axis equal vis3d;
                       
                    else % conic lens 
                        x0 = r_in( :, 1 );
                        y0 = r_in( :, 2 );
                        z0 = r_in( :, 3 );
                        e1 = e( :, 1 );
                        e2 = e( :, 2 );
                        e3 = e( :, 3 );
                        k = surf.k;
                        a = 1 + k;
                        R = surf.R;
                        
                        if a == 0 % paraboloid, special case
                            A = e2.^2 + e3.^2;
                            B = e1 * R - e2 .* y0 - e3 .* z0;
                            D = B.^2 - A .* ( -2 * R * x0 + y0.^2 + z0.^2 );
                            D( D < 0 ) = 0; %  mark no intersection as 0, I is nulled below anyway
                            
                            d01 = ( y0.^2 + z0.^2 ) / ( 2 * R ) - x0; % distance to the intersection for the ray || to the paraboloid
                            d02 = d01;
                        else
                            A = e1 .* ( a * e1.^2 + e2.^2 + e3.^2 );
                            B = e1.^2 * R - a * e1.^2 .* x0 - e1 .* e2 .* y0 - e1 .* e3 .* z0;
                            D = e1.^2 .* ( e3.^2 .* ( 2 * R * x0 - a * x0.^2 - y0.^2 ) + 2 * e2 .* e3 .* y0 .* z0 - ...
                                2 * e1 .* ( R - a * x0 ) .* ( e2 .* y0 + e3 .* z0 ) + e2.^2 .* ( 2 * R * x0 - a * x0.^2 - z0.^2 ) + ...
                                e1.^2 .* ( R^2 - a * ( y0.^2 + z0.^2 ) ) );
                            D( D < 0 ) = 0; %  mark no intersection as 0, I is nulled below anyway
                            
                            A0 = 2 * ( 2 * a * e2 .* e3 .* y0 .* z0 + e2.^2 .* ( R^2 - a * z0.^2 ) + e3.^2 .* ( R^2 - a * y0.^2 ) );
                            B0 = a * ( 1 - a ) * ( a * ( e1 .* x0 - e2 .* y0 - e3 .* z0 ) - R * e1 ) .* ( y0.^2 + z0.^2 );
                            d01 =  B0 ./ A0;   % distance to the intersection for the ray || to the hyperboloid sides
                            d02 = -B0 ./ A0;
                        end
                        
                        d1 = ( B + sqrt( D ) ) ./ A;
                        d2 = ( B - sqrt( D ) ) ./ A;
                        
                        % it is necessary to eliminate infinities before the following logical operation
                        d1(  ~isfinite( d1 )  ) = 0;
                        d2(  ~isfinite( d2 )  ) = 0;
                        d01( ~isfinite( d01 ) ) = 0;
                        d02( ~isfinite( d02 ) ) = 0;
                        
                        d( :, 1 ) = ( abs( A ) <= eps ) .* d01 + ... % ray parallel to the paraboloid / hyperboloid, special case
                            ( abs( A ) >  eps ) .* d1;
                        d( :, 2 ) = ( abs( A ) <= eps ) .* d02 + ... % ray parallel to the paraboloid / hyperboloid, special case
                            ( abs( A ) >  eps ) .* d2;
                        
                        % find the shortest positive distance to the (two) intersections along the ray
                        if a * R < 0 % hyperboloid with positive R, or ellipsoid with negative R: we want to consider only the top branch
                            ind = x0 - R / a * ( 1 + sqrt( 1 - a * ( y0.^2 + z0.^2 ) / R^2 ) ) < 0; % sources below the bottom branch
                            d( ind ) = -d( ind ); % disregard intersections with the lower branch by setting the distances to negative values
                        elseif a * R > 0 % hyperboloid with negative R, or ellipsoid with positive R: we want to consider only the bottom branch
                            ind = x0 - R / a * ( 1 - sqrt( 1 - a * ( y0.^2 + z0.^2 ) / R^2 ) ) > 0; % sources above the bottom branch
                            d( ind ) = -d( ind ); % disregard intersections with the lower branch by setting the distances to negative values
                        end
                        d( d < 0 ) = realmax; % intensities for these rays (non-intersecting the surface) will be set to 0 anyway
                        
                        [ ~, ii ] = min( d, [], 2 );
                        ind = sub2ind( size( d ), ( 1:size( d, 1 ) )', ii ); % linear index of the min value
                        d = abs( d( ind ) );
                        
                        % form the intersection vector
                        rinter = r_in + repmat( d, 1, 3 ) .* e;
                        
                        % mark intersections outside the lens diameter as 0 intensity
                        out = sum( rinter( :, 2 : 3 ).^2, 2 ) >= ( surf.D / 2 )^2;
                        rays_out.I( out, : ) = 0;
                        rinter( out, : ) = Inf; % put to infinity to prevent drawing these rays
                        
                        % find normals
                        r2yz = ( rinter( :, 2 ).^2 + rinter( :, 3 ).^2 ) / R^2; % distance to the lens center along the lens plane in units of lens R
                        if a == 0 % parabola, special case
                            c = 1 ./ sqrt( 1 + r2yz );
                            s = sqrt( 1 - c.^2 );
                        else
                            s = sqrt( r2yz ) ./ sqrt( 1 - k * r2yz );
                            c = sqrt( 1 - s.^2 );
                        end
%                         if R > 0 % add corrugations
%                             amp = 0.001;
%                             per = 1;
%                             c = 1 ./ sqrt( 1 + ( sqrt( r2yz ) ./ sqrt( 1 - ( 1 + k ) * r2yz ) + 2 * pi * amp / per * sin( 2 * pi / per * R * sqrt( r2yz ) ) ).^2 );
%                             s = sqrt( 1 - c.^2 );
%                         end
                        th = atan2( rinter( :, 3 ), rinter( :, 2 ) ); % rotation angle to bring r into XZ plane
                        en = -sign( R ) * [ -sign( R ) * c, s .* cos( th ), s .* sin( th ) ]; % make normal sign positive wrt ray
                    
                        if isa( surf, 'Retina' ) % calculate retinal image
                            % scale to a unit spherical surphace
                            rtr = rinter .* ( 1 + surf.k ) / surf.R;
                            rtr( :, 1 ) = rtr( :, 1 ) - 1;
                            [ az, el ] = cart2sph( rtr( :, 2 ), rtr( :, 3 ), rtr( :, 1 ) ); % YZX to account for Optometrika's coordinate system
                            surf.image = hist2( az, el, self.I, ...
                                linspace( -pi, pi, surf.azbins ), ...
                                linspace( -pi/2, surf.ang, surf.elbins ) );
                        end
                    end
                    
                    % return to the original RF
                    if surf.rotang ~= 0 % needs rotation
                        rays_out.r = rodrigues_rot( rinter, surf.rotax, surf.rotang );
                        nrms = rodrigues_rot( en, surf.rotax, surf.rotang );
                    else
                        rays_out.r = rinter;
                        nrms = en;
                    end
                    
                    rays_out.r = rays_out.r + repmat( surf.r, self.cnt, 1 );

                otherwise
                    error( [ 'Surface ' class( surf ) ' is not defined!' ] );
            end
        end
         
        function rays_out = interaction( self, surf )
            % INTERACTION calculates rays properties after interacting with
            % a Surface
            
            % find intersections and set outcoming rays starting points
            [ rays_out, nrms ] = self.intersection( surf );
            
            % calculate refraction
            switch( class( surf ) )
                case { 'Aperture', 'Screen', 'Retina' }
                    % nothing to be done
                    
                case { 'GeneralLens' 'FresnelLens' 'Plane' 'Lens' }
                    % calculate refraction (Snell's law)
                    med1 = surf.glass{1};
                    med2 = surf.glass{2};

                    cs1 = dot( nrms, self.n, 2 );
                    if isa( surf, 'FresnelLens' ) % remove rays coming to the Fresnel surface from the inside (through the cylindrical walls).
                        bads = cs1 < 0;
                        rays_out.I( bads ) = 0;
                        rays_out.r( bads, : ) = Inf * rays_out.r( bads, : );
                    end
                    
                    if strcmp( med1, 'mirror' ) || strcmp( med2, 'mirror' ) % if a mirror
                        rays_out.n = self.n - 2 * repmat( cs1, 1, 3 ) .* nrms; % Snell's law of reflection
                        
                        miss = rays_out.I == 0; % indices of the rays missing the mirror
                        % restore these rays
                        rays_out.I( miss ) = self.I( miss );
                        rays_out.r( miss, : ) = self.r( miss, : );
                        rays_out.n( miss, : ) = self.n( miss, : );
                        if strcmp( med1, 'mirror' ) && strcmp( med2, 'air' ) % mirror facing away
                            rays_out.I( cs1 > 0 & ~miss ) = 0; % zero rays hitting such mirror from the back
                        elseif strcmp( med1, 'air' ) && strcmp( med2, 'mirror' ) % mirror facing toward me
                            rays_out.I( cs1 < 0 & ~miss) = 0; % zero rays hitting such mirror from the back
                        end
                    else
                        if sum( self.w - self.w( 1 ) ) ~= 0 % not all wavelengths are the same
                            n_refr_in  = refrindx( self.w, med1 );
                            n_refr_out = refrindx( self.w, med2 );
                        else
                            n_refr_in  = repmat( refrindx( self.w( 1 ), med1 ), self.cnt, 1 );
                            n_refr_out = repmat( refrindx( self.w( 1 ), med2 ), self.cnt, 1 );
                        end
                        rn = n_refr_in ./ n_refr_out; % ratio of in and out refractive indices
                        cs2 = sqrt( 1 - rn.^2 .* ( 1 - cs1.^2 ) );
                        rays_out.n = repmat( rn, 1, 3 ) .* self.n - repmat( rn .* cs1 - sign( cs1 ) .* cs2, 1, 3 ) .* nrms; % refracted direction
                        
                        % calculate transmitted intensity (Fresnel formulas)
                        rs = ( rn .* cs1 - cs2 ) ./ ( rn .* cs1 + cs2 );
                        rp = ( cs1 - rn .* cs2 ) ./ ( cs1 + rn .* cs2 );
                        refraction_loss = ( abs( rs ).^2 + abs( rp ).^2 ) / 2;
                        % handle total reflection
                        tot = imag( cs2 ) ~= 0;
                        rays_out.n( tot, : ) = 0; % zero direction for such rays
                        refraction_loss( tot ) = 1;
                        out = rays_out.I == 0;
                        rays_out.I = ( 1 - refraction_loss ) .* rays_out.I; % intensity of the outcoming rays
                        % reset zero rays to zero (refraction loss can be Inf or Nan for such rays)
                        rays_out.I( out ) = 0;
                        rays_out.nrefr = n_refr_out; % refractive indices of the outcoming rays
                    end                     
                otherwise
                    error( [ 'Surface ' class( surf ) ' is not defined!' ] );
            end
        end
        
        function self = append( self, rays )
            % append rays to the current bundle
            self.r = [ self.r; rays.r ];
            self.n = [ self.n; rays.n ];
            self.w = [ self.w; rays.w ];
            self.I = [ self.I; rays.I ];
            self.nrefr = [ self.nrefr; rays.nrefr ];
            self.att = [ self.att; rays.att ];
            self.color = [ self.color; rays.color ];
            self.cnt = self.cnt + rays.cnt;                
        end
        
        function self = truncate( self )
            % remove rays with zero intensity
            ind = self.I == 0;
            self.r( ind, : ) = [];
            self.n( ind, : ) = [];
            self.w( ind, : ) = [];
            self.I( ind, : ) = [];
            self.nrefr( ind, : ) = [];
            self.att( ind, : ) = [];
            self.cnt = self.cnt - sum( ind );
        end
                            
        function [ av, dv, nrays ] = stat( self )
            % calculate mean and standard deviation of the rays startingpoints
            vis = self.I ~= 0; % visible rays
            norm = sum( self.I( vis ) );
            av = sum( repmat( self.I( vis ), 1, 3 ) .* self.r( vis, : ) ) ./ norm;
            dv = sqrt( sum( self.I( vis ) .* sum( ( self.r( vis, : ) - repmat( av, sum( vis ), 1 ) ).^2, 2 ) ) ./ norm );
            nrays = sum( vis );
        end

        function [ x0, cv, ax, ang, nrays ] = stat_ellipse( self )
            % calculate parameters of an ellips covering the projection of
            % the rays startingpoints onto YZ plane
            vis = self.I ~= 0; % visible rays
            nrays = sum( vis );
            [ x0, cv, ax, ang ] = ellipse_fit( self.r( vis, 2:3 ) );
            %[ x0, cv, ax, ang ] = ellipse_fit( self.r( vis, : ) );
        end

        function [ mu, sigma, lambda, angle, nrays ] = stat_exGaussian( self )
            % calculate parameters of the exponentially modified Gaussian distribution (EMG) fitting the rays startingpoints
            vis = self.I ~= 0; % visible rays
            nrays = sum( vis );
            [ mu, sigma, lambda, angle ] = exGaussian_fit2D( self.r( vis, 2:3 ) );
        end
        
        function [ av, dv ] = stat_sph( self, sph_pos )
            % calculate mean and standard deviation of the rays
            % startingpoints in spherical coordinates (e.g., on a retina)
            % returns average and std for [ azimuth, elevation, radius ]
            rs = self.r - sph_pos; % coordinates wrt sphere center
            sph = cart2sph( rs( :, 1 ), rs( :, 2 ), rs( :, 3 ) ); % [ az el r ]
            vis = self.I ~= 0; % visible rays
            norm = sum( self.I( vis ) );
            av = sum( repmat( self.I( vis ), 1, 3 ) .* sph( vis, 1:2 ) ) ./ norm;
            dv = sqrt( sum( self.I( vis ) .* sum( ( sph( vis, 1:2 ) - repmat( av, sum( vis ), 1 ) ).^2, 2 ) ) ./ norm );
        end
        
        function r2 = dist2rays( self, p )
            t = self.r - repmat( p, self.cnt, 1 ); % vectors from the point p to the rays origins
            proj = t * self.n'; % projections of the vectors t onto the ray
            r2 = sum( ( t - repmat( proj, 1, 3 ) .* self.n ).^2 );
        end
        
        function f = focal_point( self )
            ind = self.I ~= 0;
            sn = self.n( ind, : );
            sr = self.r( ind, : );
            scnt = sum( ind );
            nav = sum( sn ); % average bundle direction
            nav = nav / sqrt( sum( nav.^2, 2 ) ); % normalize the average direction vector
            rav = mean( sr ); % average bundle origin
            if exist( 'fminunc', 'file' ) % requires optimization toolbox
                options = optimoptions( 'fminunc', 'Algorithm', 'quasi-newton', 'Display', 'off', 'Diagnostics', 'off' );
                plr = rav + nav * fminunc( @scatter, 0, options, sr, sn, rav, nav ); % optimal plane position
            else % no optimization toolbox
                options = optimoptions( 'fminunc', 'MaxFunEvals', 2000, 'Display', 'off', 'Diagnostics', 'off' );
                plr = rav + nav * fminsearch( @scatter, 0, options, sr, sn, rav, nav );
            end
            d = dot( repmat( nav, scnt, 1 ), repmat( plr, scnt, 1 ) - sr, 2 ) ./ ...
                dot( sn, repmat( nav, scnt, 1 ), 2 );
            % calculate intersection vectors
            rinter = sr + repmat( d, 1, 3 ) .* sn; % intersection vectors
            f = mean( rinter ); % assign focus to the mean of the intersection points
        end
   end
end


function d = dist2( l, r0, e, surf )
rend = r0 + l * e; % the ray's end
d = ( rend( :, 1 ) - surf.funch( rend( :, 2 ), rend( :, 3 ), surf.funca, 0 ) ).^2;
end

function v = scatter( x, r, n, plr0, pln ) % scatter (var) of rays intersections with a plane with norm pln and origin plr
plr = plr0 + x * pln;  % the tested plane position
cnt = size( n, 1 );
d = dot( repmat( pln, cnt, 1 ), repmat( plr, cnt, 1 ) - r, 2 ) ./ ...
    dot( n, repmat( pln, cnt, 1 ), 2 );
% calculate intersection vectors
rinter = r + repmat( d, 1, 3 ) .* n; % intersection vectors
v = sum( sum( ( rinter - repmat( mean( rinter ), cnt, 1 ) ).^2, 2 ) ); % variance of the intersection points
end

