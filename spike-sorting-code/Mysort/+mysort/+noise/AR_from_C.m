function [A, gamma_y, Sigma_u_hat] = AR_from_C(C, nC, max_order)
%     build covariance sequence for lags [0,max_order+1] from block toeplitz matrix
%     
%     Here the input matrix C is assumed to exhibit a block structure with respect to the
%     multichanneled data and the temporal lag, where the temporal per channel forms the
%     Toeplitz structure, and the multichanneled part the block structure.
%     
%     The output is a sequence of one nc x nc covariance matrix per lag in the toeplitz
%     structure of C. Also the AR coefficients for the order specified
%     
    % inits
    assert(size(C,1) == size(C,2), 'C is not square');
    tf = size(C,1) / nC;
    assert(tf == round(tf), 'nc does not match tf');
    
    order = min(tf - 1, max_order + 1);    
    assert(order>1, 'order must be at least 2');
    
    gamma_y = zeros(order, nC, nC);
    GammaY0 = zeros((order - 1) * nC, (order - 1) * nC);
    
    % build gamma_y_seq from C
    for o = 0:order-1
        for c1 = 0:nC-1
            for c2 = 0:nC-1
                gamma_y(o+1, c1+1, c2+1) = C(c1*tf+o +1, c2*tf +1);
            end
        end
    end
    
    % build GammaY0 from gamma_y_seq, blockwise
    for i = 0:order - 2
        for j = i : order - 2
            GammaY0(i*nC +1:(i+1)*nC, j*nC +1:(j+1)*nC) = squeeze(gamma_y(j-i+1,:,:));
            if j > i
                GammaY0(j*nC +1:(j+1)*nC, i*nC +1:(i+1)*nC) = squeeze(gamma_y(j-i+1,:,:))';
            end
        end
    end
    
    % build calculation matrices
    gamma_y_seq = zeros(nC, nC*(order-1));
    for i=1:order-1
        gamma_y_seq(:,(i-1)*nC+1:i*nC) = squeeze(gamma_y(i+1,:,:));
    end
    %GammaY0inv = inv(GammaY0);    
    % Compute A = (A1 A2 ... Ap) matrix by Formula 3.3.6, p.83, Lütkepohl 2005
    % A = YX'/(XX')
    A = gamma_y_seq/GammaY0;
    % This corresponds to eq.31 in Neumaier&Schneider2001, left side
    % B = (inv(R11) * R12)'; but B also containes the intercept terms

    Gamma_y0 = mysort.noise.lagCfromTimeC(C, nC, 0);
    % compute sigma u hat by eq.3.2.23, p.80 and eq. 3.4.11, p.90
    % YY' - YZ'inv(ZZ')ZY' which is (with (YZ')' = ZY'):
    % GY0 - YX'inv(XX')(YX')'
    Sigma_u_hat = Gamma_y0 - A*gamma_y_seq';
    
    
    % compute sigma u hat by eq. 3.4.11, p.90
    % A*gamma_y_seq' - gamma_y_seq*A' = 0
    
    % A*GammaY0*A' IS EXACTLY  A*gamma_y_seq'
    %Sigma_u_hat = Gamma_y0 - A*gamma_y_seq' - gamma_y_seq*A' + A*GammaY0*A';
end

