spike_times_amps_chan = [80 -80  1
                         100 20  1
                         110 -20 1
                         120 20  1
                         130 -10 1];
mergeSpikesMaxDist = 10;
minGroupSize = 3;
minAmplitude = 20;
                     
[spike_times_amps_chan keepspikes] = mysort.spiketrain.findOscillationGroups(spike_times_amps_chan, ...
    mergeSpikesMaxDist, minGroupSize, minAmplitude)

spike_times_amps_chan = [ 80 -80 1
                         100 20  1
                         110 20 1
                         120 20  1
                         140 -10 1
                         280 -80 2
                         300 20  2
                         310 -20 2
                         320 20  2
                         340 -10 2];
mergeSpikesMaxDist = 10;
minGroupSize = 3;
minAmplitude = 20;
                     
[spike_times_amps_chan keepspikes] = mysort.spiketrain.findOscillationGroups(spike_times_amps_chan, ...
    mergeSpikesMaxDist, minGroupSize, minAmplitude)

spike_times_amps_chan = [ 80 -80 1
                         100 20  1
                         110 -20 1
                         120 20  1
                         140 -10 1
                          80 -80 2
                         100 20  2
                         110 -20 2
                         120 20  2
                         140 -10 2];
mergeSpikesMaxDist = 10;
minGroupSize = 4;
minAmplitude = 20;
                     
[spike_times_amps_chan keepspikes] = mysort.spiketrain.findOscillationGroups(spike_times_amps_chan, ...
    mergeSpikesMaxDist, minGroupSize, minAmplitude)