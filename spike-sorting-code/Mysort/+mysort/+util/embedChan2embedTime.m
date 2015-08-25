function y = embedChan2embedTime(x, nC)
    % Changes the representation of x from a channel embedding to a time
    % embedding. If x is a vector the orientation is ignored. If x is a
    % matrix every ROW of x is treated as one vector
    %
    % The channel embedding keeps the different channels together:
    %
    % x_c = x(a1 b1 ... K1 a2 b2 ...)
    %
    % where "a1" is the 1st sample on channel a. 
    % Time embedding means
    %
    % x_t = x(a1 a2 ... aN b1 b2 ...)
    %
    % Here an example for four channels (a to d) and 3 timelags
    %
    % time embedding (x):
    %      a1 a2 a3  b1 b2 b3  c1 c2 c3  d1 d2 d3
    % tau: |--|--|
    % nC:     |---------|---------|---------|
    %
    % channel embedding (y):
    %      a1 b1 c1 d1  a2 b2 c2 d2  a3 b3 c3 d3
    % nC:  |--|--|--|
    % tau:      |------------|------------|
    
    y = zeros(size(x));
    if size(x,1) > 1 && size(x,2) > 1
        % x is matrix 
        Tf = size(x,2)/nC;
        for tau=0:Tf-1
            for c=1:nC
                y(:, Tf*(c-1) + tau +1) = x(:, nC*tau + c);
            end
        end
    else
        % x is vector
        Tf = length(y)/nC;

        for tau=0:Tf-1
            for c=1:nC
                y(Tf*(c-1) + tau +1) = x(nC*tau + c);
            end
        end
    end

