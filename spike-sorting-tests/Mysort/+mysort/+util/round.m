function x = round(x, digit)
    if nargin == 1
        digit = 0;
    end
    
    x = x*10^digit;
    x = round(x);
    x = x/10^digit;