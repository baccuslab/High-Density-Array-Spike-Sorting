
function [classification_types mapping] = match_clustering(c1, c2)
%     n1 = length(c1);
%     n2 = length(c2);
    
    true_class_ids = unique(c1);
    class_ids = unique(c2);
    
    ncl1 = length(true_class_ids);
    ncl2 = length(class_ids);
    
    matches = zeros(ncl1, ncl2);
    for i=1:ncl1
        tc = true_class_ids(i);
        for j=1:ncl2
            c = class_ids(j);
            matches(i,j) = sum(c1==tc & c2==c);
        end
    end
    
    mapping = zeros(1, ncl1);
    for a=1:ncl1
        [i, j] = mysort.util.matrixArgMax(matches);
        mapping(i) = j;
        matches(i,:) = -1;
        matches(:,j) = -1;
    end
    def = mysort.util.defs();
    
    classification_types = def.cl * ones(size(c1));
    for i=1:ncl1
        tc = true_class_ids(i);
        j  = mapping(i);
        if j>0
            c  = class_ids(j);
            classification_types(c1==tc & c2==c) = def.tp;
        end
%         for j=1:ncl2
%             c_ = class_ids(j);
%             if c_ ~= c
%                 classification_types(c1==tc & c2==c_) = def.fp;
%             end
%         end
    end
