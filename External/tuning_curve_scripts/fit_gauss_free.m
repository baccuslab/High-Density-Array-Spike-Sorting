function [baseline,maxfiring1,maxfiring2,PO1,PO2,sigma1,sigma2,resnorm]=...
    fit_gauss_free(tuning)
% [baseline,maxfiring1,maxfiring2,PO1,PO2,sigma1,sigma2,resnorm]=...
%    fit_gauss_free(tuning)
% Fits two (circular) Gaussians with independent heights and sigmas
% PO2 is relative to PO1 (and forced to be between 90 and 270)
%
% tuning should be a matrix with tuning curves (a row for every cell)

vDirRange=0:22.5:359.9;

% Gaussian fit
PO1=zeros(1,size(tuning,1));
PO2=zeros(1,size(tuning,1));
sigma1=zeros(1,size(tuning,1));
sigma2=zeros(1,size(tuning,1));
baseline=zeros(1,size(tuning,1));
maxfiring1=zeros(1,size(tuning,1));
maxfiring2=zeros(1,size(tuning,1));
resnorm=zeros(1,size(tuning,1));
options=optimset('lsqcurvefit');
options.Display='off';
LowerBounds=[0 0 0 -360 90 0 0]; UpperBounds=[100 200 200 360 270 180 180];
for c=1:size(tuning,1)
    %inital guess:
    [m,i]=max(tuning(c,:));
    base=min(tuning(c,:));
    x0=[base (m-base) 0.8*(m-base) vDirRange(i) 180 30 30];
    [fit,ressquared]=lsqcurvefit(@gauss_free,x0,vDirRange(:)',tuning(c,:),LowerBounds,UpperBounds,options);
    baseline(c)=fit(1);
    maxfiring1(c)=fit(2);
    maxfiring2(c)=fit(3);
    PO1(c)=fit(4);
    if (PO1(c)<0)
        PO1(c) = PO1(c) + 360;
    end
    PO2(c)=fit(5);
    sigma1(c)=fit(6);
    sigma2(c)=fit(7);
    resnorm(c)=ressquared;
end

end
