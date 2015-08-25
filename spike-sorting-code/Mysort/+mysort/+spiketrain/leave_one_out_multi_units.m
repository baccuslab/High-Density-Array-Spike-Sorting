
function mua_others = leave_one_out_multi_units(spike_trains)
    % computes for all spike trains (columns) of the input, the
    % corresponding multi unit, that is given by the sum of the other spike
    % trains.
    % Input:
    %    spike_trains - cell array with trials as rows and units als
    %                   columns
    % Output: 
    %    mua_others   - cell array with trials as rows and one column per
    %                   input spike train. Ever column represents the mulit
    %                   unit spike trains of all other units than this one
    mua_others = {};
    for unit1=1:size(spike_trains,2)
        mua_others{1, unit1} = {};
        for trial=1:size(spike_trains,1)
            mua_others{trial, unit1} = [];
            for unit2 = 1:size(spike_trains,2)
                if unit1~=unit2
                    mua_others{trial, unit1} = [mua_others{trial, unit1} spike_trains{trial, unit2}'];
                end
            end
            mua_others{trial, unit1} = sort(mua_others{trial, unit1});
        end
    end