function h = upsampleDesignFilter(p, N, bta)
% Designs the filter for resampling, based on the code found in the
% standard matlab resample function

if nargin < 4,  bta = 5;  end   %--- design parameter for Kaiser window LPF
if nargin < 3,   N = 10;   end
assert(N>0, 'N must be greater than 0!');

% design filter
fc = 1/2/p;
L = 2*N*p + 1;
h = p*firls( L-1, [0 2*fc 2*fc 1], [1 1 0 0]).*kaiser(L,bta)' ;

Lhalf = (L-1)/2;
