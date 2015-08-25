function c = lagCfromTimeC(Cte, nC, lag)
    error('this function is depricated!');
    % return submatrix for specific lag from an time embedded (te) covariance
    % matrix C with nC channels
    c = zeros(nC, nC);
    for i_=1:nC
        for j_=1:nC
            c(i_,j_) = Cte((i_-1)*nC +lag+1, (j_-1)*nC +lag+1);
        end
    end
end