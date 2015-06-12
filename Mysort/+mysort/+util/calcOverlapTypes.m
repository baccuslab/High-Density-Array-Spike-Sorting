
function [overIdx overPriors] = calcOverlapTypes(overTaus, priors)
    nF = length(priors);
    overIdx = struct([]); overPriors=[];
    count = 1;
    for f1 = 1:nF
        for f2 = f1+1:nF
            for t=1:length(overTaus)
                overIdx(count).f1 = f1;
                overIdx(count).f2 = f2; 
                tau = overTaus(t);
                overIdx(count).tau = tau; 
                overIdx(count).tau_f1  = floor(tau/2); 
                overIdx(count).tau_f2  = -ceil(tau/2); 
                overPriors(count,1) = priors(f1)*priors(f2);  
                count = count +1;
            end
        end
    end   
end