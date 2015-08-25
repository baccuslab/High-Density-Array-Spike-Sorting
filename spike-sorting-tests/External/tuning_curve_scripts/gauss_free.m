function y=gauss_free(X,XData)
% gauss_free(X,XData)
% 2 Circular Gaussians with baseline, two different sigmas allowed, peaks
% do not have to be 180 degrees apart
% PD2 is relative to PD1
% X=[baseline maxAmplitude1 maxAmplitude2 PD1 PD2 stddev1 stddev2]

X=X(:)';
X=X(:,:,ones(size(XData,2),1));
XData=reshape(XData,size(XData,1),1,size(XData,2));
y=squeeze(X(:,1,:) + X(:,2,:).*exp(-(mod(180+XData-X(:,4,:),360)-180).^2./(2.*X(:,6,:).^2))+...
    X(:,3,:).*exp(-(mod(180+XData-(X(:,4,:)+X(:,5,:)),360)-180).^2./((2.*X(:,7,:).^2))));