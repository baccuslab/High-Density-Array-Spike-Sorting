function example7()
%
% test a Fresnel lens object
%
% Copyright: Yury Petrov, Oculus VR, 01/2014
%

% create a container for optical elements (Bench class)
bench = Bench;

% add optical elements in the order they are encountered by light rays

% front lens surface

% % This is the simulated aspheric lens surface
% R = 25; % lens radius of curvature
% k = -2; % lens aspheric constant
% lens1 = Lens( [ 30 0 0 ], 58, R, k, { 'air' 'pmma' } ); % concave hyperbolic surface

% % This simulates the above apsheric lens with a Fresnel lens
% rads = linspace( 0, 58 / 2, 1000 );
% a = 1 + k;
% d = sqrt( 1 - a * ( rads / R ).^2 );
% sags = R / a * ( 1 - d ); % conic profile
% angs = 90 - 180 / pi * atan( rads ./ ( R * d ) ); % cone angles in degrees with respect to the x-axis direction
% rads = rads( 2 : end ); % leave out 0 radius
% sags = sags( 1 : end - 1 ); % leave out values after the last ring
% angs = angs( 1 : end - 1 );
% lens1 = FresnelLens( [ 30 0 0 ], rads, sags, angs, { 'air' 'pmma' } ); % Fresnel surface

% 130 cones Fresnel lens surface simulating a parabolic lens surface with R = 25.32
load( 'FresnelLens_pars.mat' );
angs = angs / 180 * pi; % transform to radians
lens1 = FresnelLens( [ 40 0 0 ], rads, sags, angs, { 'air' 'pmma' } ); % Fresnel surface


% back lens surface
lens2 = Lens( [ 47.024 0 0 ], 50.8, -114.898, 0, { 'pmma' 'air' } ); % concave hyperbolic surface

bench.append( { lens1, lens2 } );

% screen
screen = Screen( [ 82.13 0 0 ], 10, 10, 128, 128 );
bench.append( screen );

% create collimated rays with some slant
nrays = 500;
rays_in = Rays( nrays, 'collimated', [ 0 0 0 ], [ 1 0 0 ], 54, 'hexagonal' );
% rays_in = Rays( 1000, 'source', [ 30 0 0 ], [ 1 1 0 ], 1, 'linear' );

tic;
fprintf( 'Tracing rays... ' );
rays_through = bench.trace( rays_in );    % repeat to get the min spread rays

% draw bench elements and draw rays as arrows
bench.draw( rays_through, 'arrows' );  % display everything, the other draw option is 'lines'

% get the screen image in high resolution
nrays = 10000;
rays_in = Rays( nrays, 'collimated', [ 0 0 0 ], [ 1 0 0 ], 54, 'hexagonal' );
rays_through = bench.trace( rays_in );
figure( 'Name', 'Image on the screen', 'NumberTitle', 'Off' );
imshow( kron( screen.image, ones( 3 ) ), [] );

toc;

% find the tightest focus position
f = rays_through( end - 1 ).focal_point;
fprintf( 'Bundle focal point: [ %.3f %.3f %.3f ]\n', f );

end
