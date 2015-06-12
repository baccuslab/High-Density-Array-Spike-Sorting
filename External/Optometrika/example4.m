function example4()
%
% test a lens with the cosine surface profile defined in coslens.m
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
lens1 = GeneralLens( [ 40 0 0 ], 58, 'coslens', { 'air' 'bk7' }, 10, 116 ); % cosine surface
% back lens surface
lens2 = Lens( [ 58 0 0 ], 58, -50, -3, { 'bk7' 'air' } ); % concave hyperbolic surface
bench.append( { lens1, lens2 } );

% screen
screen = Screen( [ 93 0 0 ], 20, 15, 128 * 20/15, 128 );
bench.append( screen );

% create collimated rays with some slant
nrays = 500;
rays_in = Rays( nrays, 'collimated', [ 0 0 0 ], [ 1 -0.1 0 ], 58, 'hexagonal' );

tic;
fprintf( 'Tracing rays... ' );
rays_through = bench.trace( rays_in );    % repeat to get the min spread rays

% draw bench elements and draw rays as arrows
bench.draw( rays_through, 'arrows' );  % display everything, the other draw option is 'lines'

% get the screen image in high resolution
nrays = 10000;
rays_in = Rays( nrays, 'collimated', [ 0 0 0 ], [ 1 -0.1 0 ], 58, 'hexagonal' );
bench.trace( rays_in );
figure( 'Name', 'Image on the screen', 'NumberTitle', 'Off' );
imshow( kron( screen.image, ones( 3 ) ), [] );

toc;
end
