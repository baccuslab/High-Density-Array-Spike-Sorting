classdef FresnelLens < Surface
    % FRESNELLENS Implements a Fresnel lens surface given by a rotation of
    % a piecewize linear curve defining a set of cone segments. These are 
    % characterized by radii, sag points, and half-angles measured with
    % respect to the lens axis of rotation.
    %
    % Member functions:
    %
    % l = FresnelLens( r, rad, sag, the, glass ) - object constructor
    % INPUT:
    % r - 1x3 position vector
    % rad - ncones x 1 vector of outer radii for the Fresnel cones
    % sag - ncones x 1 vector of sag values for the inner cone edge
    % the - ncones x 1 vector of cone half-angles, radians
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
        ncones = 0; % number of Fresnel cones
        rad = []   % cone radii
        sag = []   % cone sags
        the = []  % cone half-angles
     end
    
    methods
        function self = FresnelLens( ar, aR, as, ath, aglass )
            if nargin == 0
                return;
            end
            if size( aR, 1 ) < size( aR, 2 )
                aR = aR';
            end
            if size( as, 1 ) < size( as, 2 )
                as = as';
            end
            if size( ath, 1 ) < size( ath, 2 )
                ath = ath';
            end
            if size( aR, 1 ) ~= size( as, 1 ) || ...
               size( aR, 1 ) ~= size( ath, 1 ) || ...
               size( as, 1 ) ~= size( ath, 1 )
                error( 'Fresnel lens R, s, and th vectors must have the same length!' );
            else
                self.ncones = size( aR, 1 );
            end
            % sort all vectors based on ascending radii
            [ aR, ind ] = sort( aR );
            as = as( ind );
            ath = ath( ind );
            
            self.r = ar;
            self.rad = aR;
            self.sag = as;
            self.the = ath;
            self.glass = aglass;
        end
        
        function display( self )
            fprintf( 'Position:\t [%.3f %.3f %.3f]\n', self.r );
            fprintf( 'Orientation:\t [%.3f %.3f %.3f]\n', self.n );
            fprintf( 'Diameter:\t %.3f\n', 2 * self.rad( end ) );
            fprintf( 'Number of Fresnel rings:\t %i\n', self.ncones );
            fprintf( 'Steepest slope (rad):\t %i\n', max( abs( self.the ) ) );
            fprintf( 'Material:\t %s | %s\n', self.glass{ 1 }, self.glass{ 2 } );
        end
        
        function h = draw( self )
            % DISPLAY the Fresnel lens surface
            nang = 50;
            h = zeros( self.ncones, 1 );
            for i = 1 : self.ncones
                if i == 1
                    [ x, y, z ] = cylinder( [ 0 self.rad( i ) ] , nang );
                    z( 1, : ) = z( 1, : ) + self.sag( i );
                    z( 2, : ) = self.sag( i ) + self.rad( i ) / tan( self.the( i ) );
                else
                    [ x, y, z ] = cylinder( [ self.rad( i - 1 ) self.rad( i ) ] , nang );
                    z( 1, : ) = self.sag( i ); % add sag to the first point
                    z( 2, : ) = self.sag( i ) + ( self.rad( i ) - self.rad( i - 1 ) ) / tan( self.the( i ) );
                end
                S = [ z(:) -y(:) x(:) ]; % put the cone into the Optometrika reference frame
                
                % rotate and shift
                if self.rotang ~= 0
                    S = rodrigues_rot( S, self.rotax, self.rotang );
                end
                x(:) = S( :, 1 ) + self.r( 1 );
                y(:) = S( :, 2 ) + self.r( 2 );
                z(:) = S( :, 3 ) + self.r( 3 );
                
                c = ones( size( x, 1 ), size( x, 2 ), 3 );
                h( i ) = surf( x, y, z, c, ...
                          'EdgeColor', 'none', 'FaceLighting','phong', 'FaceColor', 'interp', 'FaceAlpha', 0.5, ...
                          'AmbientStrength', 0., 'SpecularStrength', 1 ); % grey color, shiny
%                 h( i ) = surf( x, y, z, c, ...
%                           'EdgeColor', 'white', 'FaceLighting','none', 'FaceColor', 'interp', 'FaceAlpha', 0.5, ...
%                           'AmbientStrength', 0., 'SpecularStrength', 1 ); % grey color, shiny
            end
        end 
        
        function rotate( self, rot_axis, rot_angle )
            self.rotate@Surface( rot_axis, rot_angle ); % rotate the surface members
            if abs( rot_angle ) > pi/2
                self.the = pi - self.the;
            end
        end
        
    end
    
end

