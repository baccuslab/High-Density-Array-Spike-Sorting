
function idx = toIdx(epochs)
if isempty(epochs)
    idx = [];
    return
end
len = sum((epochs(:,2) - epochs(:,1))+1);
idx = zeros(1,len);
start = 1;
for i=1:size(epochs,1)
    elen = epochs(i,2) - epochs(i,1);
    idx(start:start+elen) = epochs(i,1):epochs(i,2);
    start = start + elen+1;
end