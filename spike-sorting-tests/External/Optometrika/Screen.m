classdef Screen < Surface
    % SCREEN implements a rectangular screen surface
    %   Detailed explanation goes here
    %
    % Member functions:
    %
    % p = Screen( r, w, h, wbins, hbins ) - object constructor
    % INPUT:
    % r - 1x3 position vector
    % w - width
    % h - height
    % wbins - number of bins in the horizontal direction
    % hbins - number of bins in the vertical direction
    % OUTPUT:
    % p - screen object
    %
    % p.display() - displays the screen p information
    %
    % p.draw() - draws the screen p in the current axes
    % 
    % p.rotate( rot_axis, rot_angle ) - rotate the screen p
    % INPUT:
    %   rot_axis - 1x3 vector defining the rotation axis
    %   rot_angle - rotation angle (radians)
    %
    % Copyright: Yury Petrov, Oculus VR, 01/2014
    %
    
    properties
        h = 1   % height
        w = 1   % width
        hbins = 128 % number of bins along y-axis
        wbins = 128 % number of bins along x-axis
        image = []; % image on the screen
    end
    
    methods
        function self = Screen( ar, aw, ah, awbins, ahbins )
            if nargin == 0
                return;
            end
            self.r = ar;
            self.h = ah;
            self.w = aw;
            self.hbins = round( ahbins );
            self.wbins = round( awbins );
        end
        
        function display( self )
            % describe self
            fprintf( 'Position:\t [%.3f %.3f %.3f]\n', self.r );
            fprintf( 'Orientation:\t [%.3f %.3f %.3f]\n', self.n );
            fprintf( 'Width:\t %.3f\n',  self.w );
            fprintf( 'Height:\t %.3f\n', self.h );
            fprintf( 'Width bins:\t %i\n',  self.wbins );
            fprintf( 'Height bins:\t %i\n', self.hbins );
            fprintf( 'Rotation axis:\t [%.3f %.3f %.3f]\n', self.rotax );
            fprintf( 'Rotation angle:\t %.3f\n',  self.rotang );
        end
        
        function hndl = draw( self )
            % draw self
            y = linspace( -self.w/2, self.w/2, self.wbins );
            z = linspace( -self.h/2, self.h/2, self.hbins );
            [ y, z ] = meshgrid( y, z );
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
            if isempty( self.image )
                c = 0.2 * ones( size( z, 1 ), size( z, 2 ) );
            else
                c = self.image;
            end
            c = flipud( c ); % to get from image to matrix form
            c = fliplr( c ); % because y-axis points to the left
            hndl = surf( x, y, z, c, 'EdgeColor', 'none', 'FaceLighting','phong', 'FaceColor', 'interp', ...
                'AmbientStrength', 0., 'SpecularStrength', 0 ); % dull
            colormap summer;
        end
        
    end
    
end

