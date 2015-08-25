
% What do want to achieve? Optimally, only one detection event per spike.
% The timepoint of this detection event should be at the negative peak?
% Optimally, two detection events for overlaps. 
%
%        +              +
% c1 -------------------------------------------------------------
%          -          -            -   -
%        +               +
% c2 -------------------------------------------------------------
%          -           -                -
%        +                +
% c3 -------------------------------------------------------------
%          -            -           -    -
%
% O:       |          |            |   |


spike_times_amps_chan = [100 -20 1   % single spike on 3 channels
                         100 -20 1
                         100 -1  2
                         100 -21 3
                         101 -15 2
                         101 -10 3
                         
                         200  10  1   % axonal spike on 1 channel
                         209 -10  1];
                                       
mergeSpikesMaxDist = 10;

[spike_times_amps_chan keepspikes] = ...
    mysort.spiketrain.mergeSingleElectrodeDetectedSpikes(...
        spike_times_amps_chan, mergeSpikesMaxDist)