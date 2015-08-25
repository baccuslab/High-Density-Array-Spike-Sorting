 classdef Lens < Surface
    % LENS Implements a lens surface given by a rotation of a conic curve
    % (conic) lens surface given by
    % z = 1/R * r^2 ./ ( 1 + sqrt( 1 - ( 1 + k ) * (r/R)^2 ) ), where
    % R is the tangent sphere radius, and k is the aspheric factor:
    % 0 < k - oblate spheroid
    % k = 0 - sphere
    % -1 < k < 0 - plolate spheroid
    % k = -1 - parabola
    % k < -1 - hyperbola
    %
    % Member functions:
    %
    % l = Lens( r, D, R, k, glass ) - object constructor
    % INPUT:
    % r - 1x3 position vector
    % D - diameter
    % R - tangent sphere radius
    % k - conic coefficient
    % glass - 1 x 2 cell array of strings, e.g., { 'air' 'acrylic' }
    % OUTPUT:
    % l - lens surface object
    %
    % l.display() - displays the surface l information
    %
    % l.draw() - draws the surface l in the current axes
    %
    % l.rotate( rot_axis, rot_angle ) - rotate the surface l
    % INPUT:
    %   rot_axis - 1x3 vector defining the rotation axis
    %   rot_angle - rotation angle (radians)
    % 
    % Copyright: Yury Petrov, Oculus VR, 01/2014
    %
    
    properties
        D = 1    % lens diameter
        k = 0    % aspheric parameter (spherical surface by default)
     end
    
    methods
        function self = Lens( ar, aD, aR, ak, aglass )
            if nargin == 0
                return;
            end
            self.r = ar;
            self.D = aD;
            self.R = aR;
            self.k = ak;
            self.glass = aglass;
            if ( self.D / 2 / self.R )^2 * ( 1 + self.k ) > 1
                % error( 'Lens Diameter is too large for its radius and apsheric parameter!' );
                self.D = -1; % signal bad parameters
            end
        end
        
        function display( self )
            fprintf( 'Position:\t [%.3f %.3f %.3f]\n', self.r );
            fprintf( 'Orientation:\t [%.3f %.3f %.3f]\n', self.n );
            fprintf( 'Diameter:\t %.3f\n', self.D );
            fprintf( 'Curv. radius:\t %.3f\n', self.R );
            fprintf( 'Asphericity:\t %.3f\n', self.k );
            fprintf( 'Material:\t %s | %s\n', self.glass{ 1 }, self.glass{ 2 } );
        end
        
        function h = draw( self )
            % DISPLAY the lens surface
            nrad = 50;
            rad = linspace( 0, self.D / 2, nrad );
            nang = 50;
            ang = linspace( 0, 2 * pi, nang );
            [ ang, rad ] = meshgrid( ang, rad );
            
            [ y, z ] = pol2cart( ang, rad );
            r2yz = ( y.^2 + z.^2 ) / self.R^2;
            
            a = 1 + self.k;
            if a == 0 % paraboloid, special case
                x = r2yz * self.R / 2;
            else
                x = self.R * r2yz ./ ( 1 + sqrt( 1 - a * r2yz ) );
            end

            S = [ x(:) y(:) z(:) ];
            
            % rotate and shift
            if self.rotang ~= 0
                S = rodrigues_rot( S, self.rotax, self.rotang );
            end
            x(:) = S( :, 1 ) + self.r( 1 );
            y(:) = S( :, 2 ) + self.r( 2 );
            z(:) = S( :, 3 ) + self.r( 3 );
            
            c = ones( size( x, 1 ), size( x, 2 ), 3 );
            h = surf( x, y, z, c, ...
                'EdgeColor', 'none', 'FaceLighting','phong', 'FaceColor', 'interp', 'FaceAlpha', 0.5, ...
                'AmbientStrength', 0., 'SpecularStrength', 1 ); % grey color, shiny
        end
        
    end
    
end

