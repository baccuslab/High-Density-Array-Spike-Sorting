function [ h, hf, hb, V ] = lens_thickness( D, Rf, Rb, kf, kb, h )
% LENS_THICKNESS calculates a lens thickness given its surface parameters
%   
% INPUT:
%    D - lens diameter
%    Rf - radius of the tangent sphere for the front surface
%    Rb - radius of the tangent sphere for the back surface
%    kf - asphreicity coefficient for the front surface
%    kb - asphreicity coefficient for the back surface
%    h  - lens height at the equator (cylindrical part)
%
% OUTPUT:
%    h - total lens thickness
%    hf - front part thickness
%    hb - back part thickness
%    V  - lens volume
%
% Copyright: Yury Petrov, Oculus VR, 01/2014
%

% calculate lens thickness
if kf == -1 % paraboloid
    hf = D^2 / ( 8 * Rf );
    Vf = pi / 8 * D^2 * hf; % parabolic volume
else % spheroids, hyperboloids
    if isinf( Rf )
        hf = 0;
        Vf = 0;
    else
        x = Rf / ( 1 + kf );
        hf = sign( Rf ) * abs( x * ( 1 - sqrt( 1 - ( 1 + kf ) * D^2 / ( 4 * Rf^2 ) ) ) );
        Vf = pi * D^2 / 8 * hf * ( 1 - hf / ( 6 * x + 3 * hf ) );
    end
end

if kb == -1 % paraboloid
    hb = D^2 / ( 8 * Rb );
    Vb = pi / 8 * D^2 * hb; % parabolic volume
else % spheroids, hyperboloids
    if isinf( Rb )
        hb = 0;
        Vb = 0;
    else
        x = Rb / ( 1 + kb );
        hb = -sign( Rb ) * abs( x * ( 1 - sqrt( 1 - ( 1 + kb ) * D^2 / ( 4 * Rb^2 ) ) ) );
        Vb = pi * D^2 / 8 * hb * ( 1 - hb / ( 6 * x + 3 * hb ) );
    end
end

Ve = pi * D^2 / 4 * h; % cylindrical volume at the equator

h = hf + hb + h; % total lens thickness
V = Vf + Vb + Ve;

end

