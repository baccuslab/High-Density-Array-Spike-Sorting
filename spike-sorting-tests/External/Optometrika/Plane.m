classdef Plane < Surface
    % PLANE implements a planar refracting or reflecting surface, e.g. one
    % face of a prism or a mirror
    %
    % Member functions:
    %
    % p = Plane( r, D, glass ) - circular planar surface constructor
    % INPUT:
    % r - 1x3 position vector
    % D - surface diameter
    % glass - 1 x 3 cell array of two strings, e.g., { 'air' 'acrylic' }
    % OUTPUT:
    % p - plane object
    %
    % p = Plane( r, w, h, glass ) - rectangular planar surface constructor
    % INPUT:
    % r - 1x3 position vector
    % w - width
    % h - height
    % glass - 1 x 3 cell array of two strings, e.g., { 'air' 'acrylic' }
    % OUTPUT:
    % p - plane object
    %
    % p.display() - displays the plane p information
    %
    % p.draw() - draws the plane p in the current axes
    % 
    % p.rotate( rot_axis, rot_angle ) - rotate the plane p
    % INPUT:
    %   rot_axis - 1x3 vector defining the rotation axis
    %   rot_angle - rotation angle (radians)
    %
    % Copyright: Yury Petrov, Oculus VR, 01/2014
    %
   
    properties
        w = 0   % width of the surface
        h = 0   % height of the surface
        D = 0   % diameter of the surface
    end
    
    methods
        function self = Plane( varargin ) % some function overloading here         
            if nargin == 0
                return;
            elseif nargin == 3
                self.r = varargin{1};
                self.D = varargin{2};
                self.glass = varargin{3};
            elseif nargin == 4
                self.r = varargin{1};
                self.w = varargin{2};
                self.h = varargin{3};
                self.glass = varargin{4};
            end
        end
                
        function display( self )
            % describe self
            fprintf( 'Position:\t [%.3f %.3f %.3f]\n', self.r );
            fprintf( 'Orientation:\t [%.3f %.3f %.3f]\n', self.n );
            if ~( self.w == 0 && self.h == 0 )
                fprintf( 'Width:\t %.3f\n',  self.w );
                fprintf( 'Height:\t %.3f\n', self.h );
            elseif self.D ~= 0
                fprintf( 'Diameter:\t %.3f\n', self.D );
            end
            fprintf( 'Rotation axis:\t [%.3f %.3f %.3f]\n', self.rotax );
            fprintf( 'Rotation angle:\t %.3f\n',  self.rotang );
            fprintf( 'Material:\t %s | %s\n', self.glass{ 1 }, self.glass{ 2 } );
        end
        
        function hndl = draw( self )
            % draw self
            if ~( self.w == 0 && self.h == 0 )
                y = [-self.w/2 self.w/2 ];
                z = [-self.h/2 self.h/2 ];
                [ y, z ] = meshgrid( y, z );
            elseif self.D ~= 0
                nrad = 50;
                rad = linspace( 0, self.D / 2, nrad );
                nang = 50;
                ang = linspace( 0, 2 * pi, nang );
                [ ang, rad ] = meshgrid( ang, rad );                
                [ y, z ] = pol2cart( ang, rad );
            end
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
            c = ones( size( x, 1 ), size( x, 2 ), 3 );
            hndl = surf( x, y, z, c, 'EdgeColor', 'none', 'FaceLighting','phong', 'FaceColor', 'interp', 'FaceAlpha', 0.5, ...
                'AmbientStrength', 0, 'SpecularStrength', 1 ); % shiny gray surface
        end
        
    end
    
end

