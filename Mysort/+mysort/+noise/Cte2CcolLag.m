function CColLag = Cte2CcolLag(Cte, nC, lag)
% Extracts the channel embedded block matrix for a specific lag from the
% time embedded (te) block matrix Cte with nC channels

Tf = size(Cte, 1)/nC;
maxLag = Tf-1;
assert(lag>=0 & lag <=maxLag, 'lag out of bounds!');
CColLag = zeros(nC, nC);

for c1=1:nC
    for c2 =1:nC
        % DO add lag only on ONE of the indices !!
        CColLag(c1,c2) = Cte((c1-1)*Tf+1 +lag, (c2-1)*Tf+1);
    end
end




