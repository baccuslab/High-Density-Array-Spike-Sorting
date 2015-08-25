classdef Aperture < Surface
    % APERTURE defines a circular opening
    %   Detailed explanation goes here
    %
    % Member functions:
    %
    % a = Aperture( r, Din, Dout ) - object constructor
    % INPUT:
    % r - 1x3 position vector
    % Din - inner diameter
    % Dout - outer diameter
    % OUTPUT:
    % a - aperture object
    %
    % a.display() - displays the aperture a information
    %
    % draw() - draws the aperture a in the current axes
    % 
    % a.rotate( rot_axis, rot_angle ) - rotate the aperture
    % INPUT:
    %   rot_axis - 1x3 vector defining the rotation axis
    %   rot_angle - rotation angle (radians)
    % 
    % Copyright: Yury Petrov, Oculus VR, 01/2014
    %
    
    properties
        D = 1      % opening diameter
        Dout = 1.5 % outside diameter
    end
    
    methods
        function self = Aperture( ar, aD, aDout )
            if nargin == 0
                return;
            end
            self.r = ar;
            self.D = aD;
            self.Dout = aDout;
            if aDout < aD
                error( 'Outer radius has to be larger than the inner radius' );
            end
        end
        
        function display( self )
            fprintf( 'Position:\t [%.3f %.3f %.3f]\n', self.r );
            fprintf( 'Orientation:\t [%.3f %.3f %.3f]\n', self.n );
            fprintf( 'Diameter in:\t %.3f\n', self.D );
            fprintf( 'Diameter out:\t %.3f\n', self.Dout );
        end
        
        function h = draw( self )
            nrad = 2;
            rad = linspace( self.D / 2, self.Dout / 2, nrad );
            nang = 30;
            ang = linspace( 0, 2 * pi, nang );
            [ ang, rad ] = meshgrid( ang, rad );
            
            y = rad .* cos( ang );
            z = rad .* sin( ang );
            x = zeros( size( y ) );
            S = [ x(:) y(:) z(:) ];
            
            % rotate and shift
            if self.rotang ~= 0
                S = rodrigues_rot( S, self.rotax, self.rotang );
            end
            x(:) = S( :, 1 ) + self.r( 1 );
            y(:) = S( :, 2 ) + self.r( 2 );
            z(:) = S( :, 3 ) + self.r( 3 );
            
            % draw
            c = 0.25 * ones( size( x, 1 ), size( x, 2 ), 3 );
            h = surf( x, y, z, c, 'EdgeColor', 'none', 'FaceLighting','phong', 'FaceColor', 'interp', ...
                'AmbientStrength', 0.25, 'SpecularStrength', 0.25 );
        end
    end
    
end

