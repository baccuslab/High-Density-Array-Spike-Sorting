function r = getRadialVectors(N, d)
    % Get N vectors pointing in equal angle distances away from the origin
    % in an d-dimensional space 
    if nargin == 1
        d = 2;
    elseif d ~=2
        error('not implemented');
    end
    
    angle = 2*pi/N;
    r = zeros(N, d);
    for i=1:N
        r(i,:) = [sin(angle*(i-1)) cos(angle*(i-1))];
    end