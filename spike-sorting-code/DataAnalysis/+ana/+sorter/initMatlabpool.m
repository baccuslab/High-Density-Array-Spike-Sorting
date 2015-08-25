function initMatlabpool(nCPUS)
if nCPUS > 1
    if matlabpool('size')==0
        disp('Opening Pool.')
        matlabpool(nCPUS);
    elseif matlabpool('size') ~= nCPUS
        disp('Re-Opening Pool.')
        matlabpool close;
        matlabpool(nCPUS);
    end
end