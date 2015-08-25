classdef Bench < handle
    % Bench class implements a system of optical elements
    % A complex optical system can be stored and manipulated as a whole by
    % making it a Bench instance.
    %
    % Member functions:
    %
    % b = Bench( obj )  - constructor function
    % INPUT:
    %   obj - an optical element, cell array of elements, or another bench
    % OUTPUT:
    %   b - bench object
    %
    % b.display() - displays bench b's information
    %
    % b.draw() - draws bench b in the current axes
    % 
    % a = b.copy() - copies bench b to bench a
    %
    % b.append( a ) - appends element a to bench b
    %
    % b.prepend( a ) - prepends element a to bench b
    %
    % b.replace( ind, a ) - replaces an element with index ind on bench b with element a
    %
    % b.remove( inds ) - removes elements located at inds on bench b
    %
    % b.rotate( rot_axis, rot_angle, rot_fl ) - rotate the bench b with all its elements
    % INPUT:
    %   rot_axis - 1x3 vector defining the rotation axis
    %   rot_angle - rotation angle (radians)
    %   rot_fl - (0, default) rotation of the bench elements wrt to the
    %   global origin, (1) rotation wrt to the bench geometric center
    %   
    % b.translate( tr_vec ) - translate the bench b with all its elements
    % INPUT:
    %   tr_vec - 1x3 translation vector
    %
    % rays_through = b.trace( rays_in ) - trace rays through optical elements
    % on the bench b
    % INPUT:
    %   rays_in - incoming rays, e.g., created by the Rays() function
    % OUTPUT:
    %   rays_through - a cell array of refracted/reflected rays of the same
    %   length as the number of optical elements on the bench.
    %
    % Copyright: Yury Petrov, Oculus VR, 01/2014
    %
  
    properties
        elem = {}     % cell array of optical elements
        cnt = 0       % counter of elements in the system
    end
    
    methods
        function self = Bench( obj )
            if nargin == 0
                return;
            end
            
            if isa( obj, 'Bench' ) % if another Bench
                obj = obj.elem;    % extract elements
            end
            
            % append object(s) to the optical system
            nobj = length( obj );
            for i = 1 : nobj;
                self.cnt = self.cnt + 1;
                if nobj == 1
                    self.elem{ self.cnt } = obj;
                elseif iscell( obj )   % other benches or cell arrays of Surfaces
                    self.elem{ self.cnt } = obj{ i };
                elseif isvector( obj ) % Rays
                    self.elem{ self.cnt } = obj( i );
                end
            end
        end
         
        function display( self )
            for i = 1 : self.cnt
                obj = self.elem{ i };
                fprintf( '\n%s:\n', class( obj ) );
                obj.display;
            end
         end
        
        function draw( self, rays, draw_fl, scale, new_figure_fl )
            % draw the bench elements and rays
            if nargin < 5
                new_figure_fl = 1; % open a new figure by default
            end
            if nargin < 4
                if nargin > 1
                    scale = ones( 1, length( rays ) );
                else
                    scale = 1;
                end
            else
                if length( scale ) == 1
                    scale = repmat( scale, 1, length( rays ) ); % make all ones
                elseif length( scale ) < length( rays )
                    if size( scale, 1 ) > size( scale, 2 )
                        scale = scale';
                    end
                    scale = [ scale ones( 1, length( rays ) - length( scale ) ) ]; % append ones 
                end
            end
            if nargin < 3
                draw_fl = 'arrows';
            end
            if nargin < 2
                rays = [];
            end
            
            if new_figure_fl == 1
                fname = dbstack;  % get debugging info
                [ ~, fname ] = fname.name; % get the second (original) call function name
                figure( 'Name', [ 'OPTOMETRIKA: ' fname ], 'NumberTitle', 'Off', ...
                    'Position', [ 0 0 1024 1024 ], ...
                    'Color', 'k' );
            end
            hold on;
            for i = 1 : self.cnt
                obj = self.elem{ i };
                obj.draw();
            end
            
            if ~isempty( rays )
                if strcmp( draw_fl, 'lines' ) % draw ray bundles as lines
                    for i = 1 : length( rays ) - 1
                        vis = rays( i ).I ~= 0; % visible rays
                        [ unique_colors, ~, ic ] = unique( rays( i ).color, 'rows' );
                        for j = 1 : size( unique_colors, 1 )
                            cvis = vis & ( ic == j );
                            plot3( [ rays( i ).r( cvis, 1 )';  rays( i + 1 ).r( cvis, 1 )' ], ...
                                   [ rays( i ).r( cvis, 2 )';  rays( i + 1 ).r( cvis, 2 )' ], ...
                                   [ rays( i ).r( cvis, 3 )';  rays( i + 1 ).r( cvis, 3 )' ], '-', 'Color', unique_colors( j, : ) );
                        end
                   end
                elseif strcmp( draw_fl, 'arrows' )
                    for i = 1 : length( rays )
                        rays( i ).draw( scale( i ) );
                    end
                end
            end
            
            axis equal vis3d off;
            %grid on;
            view( -54, 54 );
            lighting phong;
%             light('Position',[-1 1 1],'Style','infinite');
%             light('Position',[-1 0 0],'Style','infinite');
%             light('Position',[ 1 -1 0],'Style','infinite');
            rotate3d on;
        end
        
        function b = copy( self )
            b = feval( class( self ) );
            b.cnt = self.cnt;
            for i = 1 : length( self.elem )
                b.elem{ i } = self.elem{ i }.copy;
            end
        end
        
        function append( self, obj )
            if isa( obj, 'Bench' ) % if another Bench
                obj = obj.elem;    % extract elements
            end           
            % append object(s) to the optical system
            nobj = length( obj );
            for i = 1 : nobj;
                self.cnt = self.cnt + 1;
                if nobj == 1
                    self.elem{ self.cnt } = obj;
                elseif iscell( obj )   % other benches or cell arrays of Surfaces
                    self.elem{ self.cnt } = obj{ i };
                elseif isvector( obj ) % Rays
                    self.elem{ self.cnt } = obj( i );
                end
            end
        end
        
        function prepend( self, obj )
            if isa( obj, 'Bench' ) % if another Bench
                obj = obj.elem;    % extract elements
            end         
            self.elem = fliplr( self.elem ); % reverse element direction temporarily
            % prepend object(s) to the optical system
            nobj = length( obj );
            for i = nobj : -1 : 1; % append in the opposite order
                self.cnt = self.cnt + 1;
                if nobj == 1
                    self.elem{ self.cnt } = obj;
                elseif iscell( obj )   % other benches or cell arrays of Surfaces
                    self.elem{ self.cnt } = obj{ i };
                elseif isvector( obj ) % Rays
                    self.elem{ self.cnt } = obj( i );
                end
            end
            self.elem = fliplr( self.elem ); % restitute the original order
        end
        
        function replace( self, ind, obj )
            self.elem{ ind } = obj;
        end
        
         function remove( self, inds )
             if self.cnt == 0
                 error( 'The bench is already empty!' );
             else
                self.elem( inds ) = [];
                self.cnt = self.cnt - length( inds );
             end
         end
         
         function rotate( self, rot_axis, rot_angle, rot_fl )
             if nargin < 4
                 rot_fl = 0;
             end
             cntr = [ 0 0 0 ];
             if rot_fl == 1 % rotate around the geometric center of the bench
                 for i = 1 : self.cnt % loop through the optic system
                     cntr = cntr + self.elem{ i }.r;
                 end
                 cntr = cntr / self.cnt;
             end
            % rotate bench elements
            for i = 1 : self.cnt % loop through the optic system
                self.elem{ i }.rotate( rot_axis, rot_angle ); % rotate normal
                self.elem{ i }.r = cntr + rodrigues_rot( self.elem{ i }.r - cntr, rot_axis, rot_angle ); % rotate position
            end
            if abs( rot_angle ) > pi/2 % reverse order in which the elements are encountered by rays
                self.elem = fliplr( self.elem );
            end
        end

        function translate( self, tr_vec )
            % translate bench elements by tr_vec
            for i = 1 : self.cnt % loop through the optic system
                self.elem{ i }.r = self.elem{ i }.r + tr_vec; % translate position
            end
        end
       
        function rays = trace( self, rays_in )
            % trace rays through the optic system, works for imaginary
            % image planes (screen in front of the lens)
            rays( 1, self.cnt + 1 ) = Rays; % allocate cnt instances of Rays
            rays( 1 ) = rays_in;
            for i = 1 : self.cnt % loop through the optic system
                rays( i + 1 ) = rays( i ).interaction( self.elem{ i } );
            end
        end
    end
    
end

