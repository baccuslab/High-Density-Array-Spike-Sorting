
function m = v2m(v, n)
    warning('This function is depricated! Use mysort.wf.* instead!');
    % converts a vector representing a multichannel (n channel) 
    % vector with channels concatinated into a matrix with the
    % single channel vectors (of dimension tf) in its rows
    % inputs:
    %    v  - vector of length n*tf
    %    n  - number of channel
    % output:
    %    y  - matrix of dimension n x tf
    assert(n>0, 'n (the number of channel) must be bigger than 0')
    
    [n1 n2] = size(v);
    if n1 == 1
        m = do_one(v,n);
    elseif n2 == 1
        m = do_one(v',n);        
    else
        % Matrix is supplied, assume single vectors are in the rows!
        m = zeros(n1*n, n2/n);
        for i=1:n1
            idx = (i-1)*n+1 : i*n;
            m(idx,:) = do_one(v(i,:),n);
        end        
    end


    function m = do_one(v,n)
        tf = length(v)/n;
        if n==1; m=v; return;  end        
        assert(int32(tf)==tf, 'length of v must be multiple of n!')
        m = reshape(v, tf, n)';
    end     
end