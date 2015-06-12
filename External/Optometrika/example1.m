function focal = example1()
%
% test the Optometrika library 
%
%
% Copyright: Yury Petrov, Oculus VR, 01/2014
%

% create a container for optical elements (Bench class)
bench = Bench;

% add optical elements in the order they are encountered by light rays

% aperture
aper = Aperture( [ 40 0 0 ], 55, 80 );
bench.append( aper );

% front lens surface
lens1 = Lens( [ 40 0 0 ], 58, 40, -1, { 'air' 'bk7' } ); % parabolic surface
% back lens surface
lens2 = Lens( [ 60 0 0 ], 58, -70, -3, { 'bk7' 'air' } ); % concave hyperbolic surface
bench.append( { lens1, lens2 } );

% prism front surface
prism1 = Plane( [ 75 0 0 ], 40, 30, { 'air' 'acrylic' } );
prism1.rotate( [ 0 0 1 ], 0.3 ); % rotate the front prism plane by .3 radians about z-axis
% prism back surface
prism2 = Plane( [ 85 0 0 ], 40, 30, { 'acrylic' 'air' } );
prism2.rotate( [ 0 1 0 ], 0.2 ); % rotate the back prism plane by .2 radians about y-axis
bench.append( { prism1, prism2 } );

% screen
screen = Screen( [ 110 0 0 ], 20, 15, 128 * 20/15, 128 );
bench.append( screen );

%bench.rotate( [ 0 0 1 ], 0.15 );

% create some rays
nrays = 500;
rays_in = Rays( nrays, 'collimated', [ 0 0 0 ], [ 1 -0.1 0 ], 58, 'hexagonal' );

tic;

npos = 50;
dv = zeros( npos, 1 );
scr_x = linspace( lens2.r(1) + 30, lens2.r(1) + 60, npos );
for i = 1 : npos % loop over different screen distances
    screen.r(1) = scr_x( i );   % note that one can change element parameters after it was added to a Bench    
    rays_through = bench.trace( rays_in ); % trace rays    
    [ ~, dv( i ) ] = rays_through( end ).stat; % get stats on the last ray bundle
end

[ mdv, mi ] = min( dv );
focal = scr_x( mi ) - lens2.r(1);
fprintf( 'Focal distance: %.3f\n', focal );
fprintf( 'Bundle std: %.3f\n', mdv );

%find the tightest focus position
f = rays_through( end - 1 ).focal_point;
fprintf( 'Bundle focal point: [ %.3f %.3f %.3f ]\n', f );

% display focusing
figure( 'Name', 'Optical system focusing', 'NumberTitle', 'Off' );
plot( scr_x - lens2.r(1), dv, '-*' );
xlabel( 'Screen distance from the back lens surface', 'FontSize', 16 );
ylabel( 'Bundle focus (standard deviation)', 'FontSize', 16 );

% draw rays for the tightest focus
screen.r(1) = scr_x( mi );           % set distance for which the spread was minimal
rays_through = bench.trace( rays_in );    % repeat to get the min spread rays

% draw bench elements and draw rays as arrows
bench.draw( rays_through, 'arrows' );  % display everything, the other draw option is 'lines'
scatter3( f(1), f(2), f(3), 'w*' ); % draw the focal point as a white *

% get the screen image in high resolution
nrays = 10000;
rays_in = Rays( nrays, 'collimated', [ 0 0 0 ], [ 1 -0.1 0 ], 58, 'hexagonal' );
bench.trace( rays_in );
figure( 'Name', 'Image on the screen', 'NumberTitle', 'Off' );
imshow( kron( screen.image, ones( 3 ) ), [] );

toc;
end

