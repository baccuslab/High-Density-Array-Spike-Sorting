function [spike_times_amps_chan keepspikes] = findOscillationGroups(spike_times_amps_chan, ...
    mergeSpikesMaxDist, minGroupSize, minAmplitude)
    % finds temporally local groups of peaks with alternating sign and
    % a certain minimal amplitude of a minimal size and removes those from
    % spike_times_amps_chan
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
    
    if isempty(spike_times_amps_chan)
        keepspikes = [];
        return
    end
    
    % make sure that the spike times are in increasing order
    [spike_times_amps_chan orig_sorting]= sortrows(spike_times_amps_chan, 1);
    
    
%     if size(spike_times_amps_chan, 2)>2
%         uChans = unique(spike_times_amps_chan(:,3));
%         if length(uChans) > 1
%             thisSTA = cell(length(uChans), 1);
%             thisKS = cell(length(uChans), 1);
%             for i=1:length(uChans)
%                 thisSTA{i} = spike_times_amps_chan(spike_times_amps_chan(:,3)==uChans(i),:);
%                 [thisSTA{i} thisKS{i}] = mysort.spiketrain.findOscillationGroups(thisSTA{i}, ...
%                     mergeSpikesMaxDist, minGroupSize, minAmplitude);
%             end
%             spike_times_amps_chan = sortrows(cell2mat(thisSTA));
%             keepspikes = cell2mat(thisKS);         
%             return
%         end
%     end
    
    nS = size(spike_times_amps_chan,1);
    keepspikes = ones(nS, 1);
    
    % Detect local groups of detections that are close enough
    % i is pointer on list of first spike of a group
    % k is pointer on list of current spike of a group
    i = 1;
    while i <= nS
        if abs(spike_times_amps_chan(i,2)) < minAmplitude
            i=i+1;
            continue
        end
        % Make new group, start at i
        k = i;
        while k<nS && spike_times_amps_chan(k+1,1) <= spike_times_amps_chan(k,1) + mergeSpikesMaxDist
            % accept new spike into group if it is close enough
            k = k+1;
        end
        % local group found within maxdist of first spike in group. 
        % check if group has desired number of peaks of minSize to treat it
        % as artifact
        nMembers = k-i+1;
        if nMembers >= minGroupSize
            amps = abs(spike_times_amps_chan(i:k,2));
            if size(spike_times_amps_chan,2)>2
                % check every channel if it has enough peaks individually
                chans = spike_times_amps_chan(i:k,3);
                uC = unique(chans);                
                for ci = 1:length(uC)
                    if sum(amps(chans==uC(ci)) >= minAmplitude) >= minGroupSize
                        keepspikes(i:k) = 0;
                        break
                    end
                end
            elseif sum(amps>=minAmplitude) >= minGroupSize
                keepspikes(i:k) = 0;
            end            
        end
        i = k+1;
    end
    spike_times_amps_chan = spike_times_amps_chan(keepspikes==1, :);
    keepspikes(orig_sorting) = keepspikes;