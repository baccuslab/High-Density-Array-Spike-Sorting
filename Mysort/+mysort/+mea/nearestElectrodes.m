
function [I R] = nearestElectrodes(x,y,idx,k)
    % finds k nearest electrodes for electrode idx
    
    D = zeros(length(x),1);
    p = [x(idx) y(idx)];
    for i=1:length(x)
        p2 = [x(i) y(i)];
        D(i) = norm(p-p2);
    end
    R = [(1:length(x))' D];
    R = sortrows(R,2);
    R = R(2:end,:); % closest electrode is itself
    I = R(1:k, 1);