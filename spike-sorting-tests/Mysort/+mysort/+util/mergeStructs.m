function P_base = mergeStructs(P_base, P_additional)

names = fieldnames(P_additional);
for i=1:length(names)
    if isstruct( P_additional.(names{i}))
         if isfield( P_base, names{i})
             P_base.(names{i}) = mysort.util.mergeStructs(P_base.(names{i}), P_additional.(names{i}));
         else
             P_base.(names{i}) = P_additional.(names{i});
         end
    else
         P_base.(names{i}) = P_additional.(names{i});
    end
end

end