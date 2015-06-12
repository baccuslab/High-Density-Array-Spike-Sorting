function [baseline,maxf1,maxf2,PD1,PD2,k1,k2,resnorm]=fit_vonmises_free(tuning)
% [baseline,maxf1,maxf2,PD1,PD2,k1,k2,resnorm]=fit_vonmises_free(tuning)
% Fits two von Mises functions
%
% tuning should be a matrix with tuning curves (a row for every cell)

vDirRange=0:22.5:359.9;

% vonMises fit
baseline=zeros(1,size(tuning,1));
maxf1=zeros(1,size(tuning,1));
maxf2=zeros(1,size(tuning,1));
PD1=zeros(1,size(tuning,1));
PD2=zeros(1,size(tuning,1));
k1=zeros(1,size(tuning,1));
k2=zeros(1,size(tuning,1));
resnorm=zeros(1,size(tuning,1));
options=optimset('lsqcurvefit');
options.Display='off';
LowerBounds=[0 0 0 -360 90 0.1 0.1]; UpperBounds=[100 200 200 360 270 Inf Inf];
for c=1:size(tuning,1)
    %inital guess:
    [m,i]=max(tuning(c,:));
    minimum=min(tuning(c,:));
    x0=[minimum (m-minimum) 0.8*(m-minimum) vDirRange(i) 180 2 2];
    [fit,ressquared]=lsqcurvefit(@vonmises_free,x0,vDirRange(:)',tuning(c,:),LowerBounds,UpperBounds,options);
    baseline(c)=fit(1);
    maxf1(c)=fit(2);
    maxf2(c)=fit(3);
    PD1(c)=fit(4);
    if (PD1(c)<0)
        PD1(c) = PD1(c) + 360;
    end
    PD2(c)=fit(5);
    k1(c)=fit(6);
    k2(c)=fit(7);   
    resnorm(c)=ressquared;
end

end
