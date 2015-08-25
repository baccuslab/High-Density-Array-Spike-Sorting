function s = sinc0(x_)
    s = sin(pi*x_)./(pi*x_);
    s(x_==0) = 1;