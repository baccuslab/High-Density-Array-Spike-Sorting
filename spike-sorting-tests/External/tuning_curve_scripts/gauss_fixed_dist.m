function y=gauss_fixed_dist(X,XData)
% gauss_fixed_dist(X,XData)
% 2 Circular Gaussians with baseline, peaks forced to be 180° apart,
% two different sigmas allowed
% X=[baseline maxAmplitude1 maxAmplitude2 PO stddev1 stddev2]

X=X(:)';
X=X(:,:,ones(size(XData,2),1));
XData=reshape(XData,size(XData,1),1,size(XData,2));
y=squeeze(X(:,1,:) + X(:,2,:).*exp(-(mod(180+XData-X(:,4,:),360)-180).^2./(2*X(:,5,:).^2))+...
    X(:,3,:).*exp(-(mod(XData-X(:,4,:)+360,360)-180).^2./((2*X(:,6,:).^2))));