function [baseline,maxf1,maxf2,PD,k1,k2,resnorm]=fit_vonmises_fixed_dist(tuning, vDirRange)
% [baseline,maxf1,maxf2,PD,k1,k2,resnorm]=fit_vonmises_fixed_dist(tuning)
% Fits two von Mises functions with peaks 180° apart
%
% tuning should be a matrix with tuning curves (a row for every cell)

if nargin == 1
    vDirRange=0:22.5:359.9;
end
assert(size(tuning, 2) == length(vDirRange), 'the x values must match the tuning y values!');

% vonMises fit
baseline=zeros(1,size(tuning,1));
maxf1=zeros(1,size(tuning,1));
maxf2=zeros(1,size(tuning,1));
PD=zeros(1,size(tuning,1));
k1=zeros(1,size(tuning,1));
k2=zeros(1,size(tuning,1));
resnorm=zeros(1,size(tuning,1));
options=optimset('lsqcurvefit');
options.Display='off';
LowerBounds=[0 0 0 -360 0.1 0.1]; UpperBounds=[100 200 200 360 Inf Inf];
for c=1:size(tuning,1)
    %inital guess:
    [m,i]=max(tuning(c,:));
    minimum=min(tuning(c,:));
    x0=[minimum (m-minimum) 0.8*(m-minimum) vDirRange(i) 2 2];
    [fit,ressquared]=lsqcurvefit(@vonmises_fixed_dist,x0,vDirRange(:)',tuning(c,:)',LowerBounds,UpperBounds,options);
    baseline(c)=fit(1);
    maxf1(c)=fit(2);
    maxf2(c)=fit(3);
    PD(c)=fit(4);    
    if (PD(c)<0)
        PD(c) = PD(c) + 360;
    end
    k1(c)=fit(5);
    k2(c)=fit(6);
    resnorm(c)=ressquared;
end

end
