
function C = mcDataCovariance(X,Tf, mode)
    % N -> number of parallel time series (channels)
    % T -> number of data points
    [nC,T]=size(X);

    % initialize noise Cov
    C=zeros(nC*Tf,nC*Tf);
    maxLag = Tf-1;
    %normlag = T+[-maxLag:0 -1:-1:-maxLag];
    for i=1:nC
        s11 = (i-1)*Tf+1;
        s12 = s11+Tf-1;
        % Compute the block diagonal elements
        xx=xcorr(X(i,:),maxLag, 'none');
        %xx = xx./normlag;
        h1=xx(Tf:2*Tf-1);
        C(s11:s12,s11:s12)=toeplitz(h1);
        
        % Compute the off diagonal blocks
        for j=i+1:nC
            s21 = (j-1)*Tf+1;
            s22 = s21 + Tf-1;

            xx=xcorr(X(i,:),X(j,:),maxLag, 'none');
            % normalize xcorr. Take care, that the crosscorrelation 
            % is calculated with zeropadding. this means for larger lags
            % we have slightly less data to compute the crosscorrelation
            % No, dont do this, this can cause negative eigenvalues.
            %xx = xx./normlag;
            h1=xx(Tf:2*Tf-1);
            h2=xx(Tf:-1:1);
            % h1 must not be equal to h2, xcorrelation is not symmetric!
            C(s11:s12,s21:s22)=toeplitz(h1,h2);
            C(s21:s22,s11:s12)=toeplitz(h2,h1);
        end
    end
    C = C/T;
end
