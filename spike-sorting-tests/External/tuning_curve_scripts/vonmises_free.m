function y=vonmises_free(X,XData)
% vonmises_free(X,XData)
% 2 von Mises functions with baseline, different widths allowed, peaks do
% not have to be 180° apart
% PD2 is relative to PD1
% X=[baseline,maxf1,maxf2,PD1,PD2,k1,k2]

X=X(:)';
X=X(:,:,ones(size(XData,2),1));
XData=reshape(XData,size(XData,1),1,size(XData,2));
base=X(:,1,:);
A1=X(:,2,:);  A2=X(:,3,:);
PD1=X(:,4,:); PD2=X(:,4,:)+X(:,5,:);
k1=X(:,6,:);  k2=X(:,7,:);

y=squeeze(base+A1.*exp(k1.*(cosd(XData-PD1)-1))+A2.*exp(k2.*(cosd(XData-PD2)-1)));