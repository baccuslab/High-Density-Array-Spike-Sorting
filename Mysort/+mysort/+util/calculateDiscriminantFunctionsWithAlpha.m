
function [D alphas] = calculateDiscriminantFunctionsWithAlpha(Y, M, priors)
% If Y are the filter outputs of the bayes optimal (LDA) template matchers
% and M the respective cross correlation functions (confusion tensor), then
% the function computes the bayes optimal discriminant functions
% (boundaries of LDA)
    nF = size(Y,1);
    nS = size(Y,2);
    
%     if ndims(M) == 3
%         % M is the confusion tensor, we need only the zero lag values
        Tf = floor(size(M,1)/2)+1;        
%         M = diag(squeeze(M(Tf,:,:)));
%     else
%         M = mysort.util.toColVec(M);
%     end
    
    D = zeros(nF, nS);
    alphas = zeros(nF, nS);
    for f=1:nF
        alphas(f,:) = Y(f,:)/M(Tf,f,f);
        mask = alphas(f,:)<0.1;
        D(f,:) = .5*(Y(f,:).^2) ./M(Tf,f,f) + log(priors(f));
        D(f,mask) = -inf;
    end        
end   