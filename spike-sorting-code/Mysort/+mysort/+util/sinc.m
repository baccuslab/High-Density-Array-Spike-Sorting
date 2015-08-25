
function s = sinc(x_, derivative)
%     if ~exist('derivative', 'var')
%         derivative = 0;
%     end
%     
%     assert(derivative>=0 && derivative<4, 'derivative must be a positive integer smaller than 4!');
    
%     if derivative == 0
%         s = sinc(x_);
%         return
%     end
    
    if nargin == 1
        derivative = 0;
    end
    if derivative == 0
        s = mysort.util.sinc0(x_);
        return
    end
    
    s = zeros(size(x_));
    z  = x_==0;
    nz = x_~=0;
    x = x_(nz);
    
    if derivative == 1
        s(z)  = 0;
        s(nz) = cos(pi*x)./x - ...
                sin(pi*x)./pi./x.^2;
    elseif derivative == 2
        s(z)  = -pi^2/3;
        s(nz) = -pi*sin(pi*x)./x + ...
                2*sin(pi*x)./(pi*x^3) - ...
                2*cos(pi*x)/x.^2;
    elseif derivative == 3
        s(z)  = 0;
        s(nz) = 3*pi*sin(p*ix)./x.^2 - ...
                6*sin(p*ix)./(p*ix^4) - ...
                pi^2*cos(p*ix)./x + ...
                6*cos(p*ix)./x^3;
    end
    
    
    