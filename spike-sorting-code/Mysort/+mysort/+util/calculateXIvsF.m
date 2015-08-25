
function XIvsF = calculateXIvsF(T, F, nC, zeroOutFilterArtefacts)
    % XIvsF(Tf+tau, i, k) = (T(k)*Stau*F(i))
    if ndims(T) == 3
        T = mysort.wf.t2m(T);
    end
    if ndims(F) == 3
        F = mysort.wf.t2m(F);
    end
    
    Tf = size(T, 2)/nC;
    nF = size(F, 1);
    nT = size(T, 1);
    if nargin == 3; zeroOutFilterArtefacts=0; end

    XIvsF = zeros(Tf*2-1, nF, nT);
    for f=1:nF
        for x=1:nT
            XIvsF(:,f,x) = mysort.util.mcfilt(mysort.wf.v2m(T(x,:),nC),...
                                      mysort.wf.v2m(F(f,:),nC),'full');
            if zeroOutFilterArtefacts
                XIvsF(1:floor(Tf/2),f,x) = 0;
                XIvsF(end-floor(Tf/2)+1:end,f,x) = 0;
            end
        end
    end
end  