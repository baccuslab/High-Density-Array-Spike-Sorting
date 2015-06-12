function [baseline,maxfiring1,maxfiring2,PO,sigma1,sigma2,resnorm]=...
    fit_gauss_fixed_dist(tuning)
% [baseline,maxfiring1,maxfiring2,PO,sigma1,sigma2,resnorm]=...
%    fit_gauss_fixed_dist(tuning)
% Fits two (circular) Gaussians forced to be 180° apart.
%
% tuning should be a matrix with tuning curves (a row for every cell)

vDirRange=0:22.5:359.9;

% Gaussian fit
PO=zeros(1,size(tuning,1));
sigma1=zeros(1,size(tuning,1));
sigma2=zeros(1,size(tuning,1));
baseline=zeros(1,size(tuning,1));
maxfiring1=zeros(1,size(tuning,1));
maxfiring2=zeros(1,size(tuning,1));
resnorm=zeros(1,size(tuning,1));
options=optimset('lsqcurvefit');
options.Display='off';
LowerBounds=[0 0 0 -360 0 0]; UpperBounds=[100 200 200 360 180 180];
for c=1:size(tuning,1)
    %inital guess:
    [m,i]=max(tuning(c,:));
    x0=[min(tuning(c,:)) m 0.8*m vDirRange(i) 30 30];
    [fit,ressquared]=lsqcurvefit(@gauss_fixed_dist,x0,vDirRange(:)',tuning(c,:),LowerBounds,UpperBounds,options);
    baseline(c)=fit(1);
    maxfiring1(c)=fit(2);
    maxfiring2(c)=fit(3);
    PO(c)=fit(4);
    if (PO(c)<0)
        PO(c) = PO(c)+360;
    end
    sigma1(c)=fit(5);
    sigma2(c)=fit(6);
    resnorm(c)=ressquared;
end

end
