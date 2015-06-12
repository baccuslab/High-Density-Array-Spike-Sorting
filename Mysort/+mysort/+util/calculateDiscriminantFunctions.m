
function [D shifts] = calculateDiscriminantFunctions(Y, M, priors)
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
    shifts = zeros(nF,1);
    for f=1:nF
        shifts(f) =     - .5* M(Tf,f,f) + log(priors(f));
        D(f,:) = Y(f,:) - .5* M(Tf,f,f) + log(priors(f));
    end        
end   