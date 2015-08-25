function DATA = detectSpikes(DATA, INTER, CONFIG)
    fprintf('Need to redetect spikes...\n');
    DATA.bNeedMarkUselessElecetrodes = 1;
    [t a] = hdmeagui.detectSpikes(DATA, CONFIG);
    DATA.singleChannelSpikeTimes = t;
    DATA.singleChannelSpikeAmplitudes = a;
    fprintf('Done.\n');
