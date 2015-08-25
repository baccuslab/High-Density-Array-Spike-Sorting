
T = [];
for subi=1:4
    for t=1:3
        for c=1:2
            T(:,c,t,subi) = (1*subi + 10*t + 100*c)*ones(10,1);
        end
    end
end
[Tf nC nS nSubSamples] = size(T);
vT = mysort.wf.t2v(T);
size(vT)
squeeze(vT(:,:,1))
T = mysort.wf.v2t(vT, nC);
vT = mysort.wf.t2v(T);
size(vT)
squeeze(vT(:,:,1))
squeeze(vT(:,:,2))