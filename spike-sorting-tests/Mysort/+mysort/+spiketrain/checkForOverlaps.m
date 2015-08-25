
function [O nO] = checkForOverlaps(St, maxOverlapDist)
    % St - cell array of spike trains
    % initialise all spikes as non overlaps
    O =  {};
    for i=1:length(St)
        O{i} = zeros(1,length(St{i}));
    end
    nO = zeros(1,length(St));
    if length(St) <= 1
        % only one spike train in St, thus no overlaps
        return
    end

    for i=1:length(St)
        for j=i+1:length(St)
            idxI = 1;
            idxJ = 1;
            while(idxI <= length(St{i}) && idxJ <= length(St{j}))
                if abs(St{i}(idxI) - St{j}(idxJ)) <= maxOverlapDist
                    % These two spikes are overlaps!
                    % Check if this was not already detected
                    if O{i}(idxI) == 0
                        O{i}(idxI) = 1;
                        nO(i) = nO(i) +1;
                    end

                    if O{j}(idxJ) == 0
                        O{j}(idxJ) = 1;
                        nO(j) = nO(j) +1;
                    end                        
                end

                if St{i}(idxI) < St{j}(idxJ)
                    idxI = idxI +1;
                else
                    idxJ = idxJ +1;
                end
            end
        end
    end
end