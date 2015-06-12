
function [tp, fp, fn, alignment, classifications1,classifications2] = ...
                    compareTwo(st1, st2, maxJitter)
tp=0; fp=0; fn=0; alignment = [];
classifications1 =zeros(length(st1),1);
classifications2 =zeros(length(st2),1);
def = mysort.util.defs();
i = 1;
j = 1;
while i<=length(st1) && j<=length(st2)
    if st1(i) <= st2(j)+maxJitter && ...
       st1(i) >= st2(j)-maxJitter
        tp = tp + 1;
        alignment(tp,:) = [i;j];
        classifications1(i) = def.tp;
        classifications2(j) = def.tp;  
        i=i+1;
        j=j+1;        
    elseif st1(i) < st2(j)
        fn = fn+1;
        classifications1(i) = def.fn;            
        i=i+1;        
    else
        fp = fp+1;
        classifications2(j) = def.fp;      
        j=j+1;        
    end
end
