function [spike_times_amps_chan keepspikes] = mergeSingleElectrodeDetectedSpikes(spike_times_amps_chan, mergeSpikesMaxDist)
    % merges spike trains that come from a spike detection on every
    % individual electrode of a multielectrode recording across all electrodes.
    %
    % input:
    %   spike_times_amps_chan - matrix containing one row per detection
    %                           event. The first row contains the detection
    %                           times in samples, the second row contains
    %                           the amplitude and the third row the channel
    %                           on which the detection occurred. 
    %                           The third row is optional.
    %   mergeSpikesMaxDist    - distance in samples between two detection
    %                           event to be merged
    % output:
    %   spike_times_amps_chan - spike_times_amps_chan same as the input,
    %                           just with the merged events
    %   keepspikes            - idx into the input spike_times_amps_chan
    %                           which events were kept
    
    % make sure that the spike times are in increasing order
    [spike_times_amps_chan orig_sorting]= sortrows(spike_times_amps_chan, 1);
    
    nS = size(spike_times_amps_chan,1);
    keepspikes = zeros(nS, 1);
    
    % Detect local groups of detections that are close enough
    % i is pointer on list of first spike of a group
    % k is pointer on list of current spike of a group
    % m is pointer on current maximum of group
    i = 1;
    while i <= nS
        k = i;
        m = i;
        %while k<nS && spike_times_amps_chan(k+1,1) < spike_times_amps_chan(i,1) + mergeSpikesMaxDist
         % do not use i, this will make problems if a random peak is
         % detected before an actual spike. It will break the real spike in
         % the middle        
        while k<nS && spike_times_amps_chan(k+1,1) < spike_times_amps_chan(m,1) + mergeSpikesMaxDist
        %while k<nS && spike_times_amps_chan(k+1,1) < spike_times_amps_chan(k,1) + mergeSpikesMaxDist
            % accept new spike into group
            k = k+1;
            % check if maximums pointer has to be moved
            if spike_times_amps_chan(m,2) > 0 && spike_times_amps_chan(k,2) <0
                % yes, the old one pointed to a maximum, but we found a
                % minumum !
                m = k;
            elseif abs(spike_times_amps_chan(m,2)) <= abs(spike_times_amps_chan(k,2))
                % yes, the new peak has a larger amplitude!
                m = k;
            end
        end
        % local group found within maxdist of first spike in group. find
        % maximal negative peak (if any)
        if k==i
            % only one peak, keep it
            keepspikes(k) = 1;
        elseif any(spike_times_amps_chan(i:k,2) < 0)
            [max_neg max_neg_idx] = min(spike_times_amps_chan(i:k,2));
            keepspikes(i-1+max_neg_idx) = 1;
        else
            % no negative peak, take maximum 
            [max_plus max_plus_idx] = max(spike_times_amps_chan(i:k,2));
            keepspikes(i-1+max_plus_idx) = 1;
        end
        i = k+1;
    end
    spike_times_amps_chan = spike_times_amps_chan(keepspikes==1, :);
    keepspikes(orig_sorting) = keepspikes;