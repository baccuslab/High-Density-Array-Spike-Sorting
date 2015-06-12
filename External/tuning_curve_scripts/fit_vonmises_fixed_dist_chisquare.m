function [baseline,maxf1,maxf2,PD,k1,k2,resnorm]=fit_vonmises_fixed_dist_chisquare(tuning,tuningerror,N)
% [baseline,maxf1,maxf2,PD,k1,k2,resnorm]=fit_vonmises_fixed_dist(tuning)
% Fits two von Mises functions with peaks 180° apart, takes the variance of
% the data points into account by minimizing the chi-squared error
%
% tuning should be a matrix with tuning curves (a row for every cell)
% tuningerror should give the SEM for every datapoint
% N is the number of trials for every datapoint

vDirRange=0:22.5:359.9;

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
LowerBounds=[0 0 0 0 0.1 0.1]; UpperBounds=[100 200 200 360 Inf Inf];
for c=1:size(tuning,1)
    %inital guess:
    [m,i]=max(tuning(c,:));
    minimum=min(tuning(c,:));
    x0=[minimum (m-minimum) 0.8*(m-minimum) vDirRange(i) 2 2];
    sigma=tuningerror(c,:).*sqrt(N(c,:));
    sigma(sigma==0) = max(sigma)*0.1; % prevent zeros
    [fit,ressquared]=lsqnonlin(@(x) chi_square(x,tuning(c,:),sigma),...
        x0,LowerBounds,UpperBounds,options);
    baseline(c)=fit(1);
    maxf1(c)=fit(2);
    maxf2(c)=fit(3);
    PD(c)=fit(4);    
    k1(c)=fit(5);
    k2(c)=fit(6);
    resnorm(c)=ressquared;
end

    function y=chi_square(params,data,sigma)
        %params are the parameters, the algorithm is trying to optimize
        estimation=vonmises_fixed_dist(params,vDirRange(:)');
        
        %y is the "chi" (not squared yet, lsqnonlin will minimize the sum
        %of squares of the values in y)
        y = (estimation-data)./sigma;
    end

end
