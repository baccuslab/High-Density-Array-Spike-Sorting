
tTMS = mysort.wf.v2t(S.templatesMS, S.nC);
Tf = size(tTMS,1);
upfac = 3;
tTMSup = mysort.util.resampleTensor(tTMS, upfac,1);
Tfup = size(tTMSup,1);

figure;
subplot(3,1,1)
plot(1:Tf, squeeze(tTMS(:,6,:)))
hold on
plot((1+2:Tfup+2)/upfac, squeeze(tTMSup(:,6,:)))

t1 = squeeze(tTMS(:,6,10));
t2 = squeeze(tTMS(:,6,11));

t1up = squeeze(tTMSup(:,6,10));
t2up = squeeze(tTMSup(:,6,11));

xc = sum(t1.^2) - 2*xcorr(t1, t2) + sum(t2.^2); 
xcup = sum(t1up.^2) - 2*xcorr(t1up, t2up) + sum(t2up.^2); 

xc = xcov(t1, t2, 'biased') / (std(t1)*std(t2)); 
xcup = xcov(t1up, t2up, 'biased') / (std(t1up)*std(t2up)); 


subplot(3,1,2);
plot(-(Tf-1):(Tf-1), xc); 
hold on
plot((-(Tfup-1):(Tfup-1))/3, xcup, 'g'); 


[xc ij maxis shifts] = mysort.wf.tXCorr(tTMS);
[xcup ij maxis shifts] = mysort.wf.tXCorr(tTMSup);

for i=1:size(tTMS,3)
    
    myidx = ij(:,1)==i;
    myjidx = find(myidx & maxis'>.95);
    myj = unique(ij(myjidx,2));
    if ~isempty(myj)
        mysort.plot.waveforms(tTMS(:,:,[i; myj]));
        title([num2str(i) ' - ' num2str(myj')]);
    end
end

