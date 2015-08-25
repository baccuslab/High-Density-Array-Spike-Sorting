% 'red' / Crystal Cove lens parameters:
D = 40.91;  % lens diameter
h = 5.26;     % lens thickness
glass = 'polycarbonate';  % lens material 
Rf = 24.2189;   % screen surface radius of curvature
Rb = -212.774; % eye surface radius of curvature 
kf = -0.863502; % front conic coefficient
kb = -60.6277;  % back conic coefficient
dscr = 34.9 + 8.86 + h/2; % distance to the screen from the lens equator
 
bench = Bench;
% screen (front) lens surface
lensF = Lens( [ dscr 0 0 ], D, Rf, kf, { 'air' glass } );
bench.append( lensF );
lth = lens_thickness( D, Rf, Rb, kf, kb, h );
% eye (back) lens surface
lensB = Lens( [ dscr + lth 0 0 ], D, Rb, kb, { glass 'air' } );
bench.append( lensB );
bench.draw;

