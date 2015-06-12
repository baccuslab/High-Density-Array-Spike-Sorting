st1 = [ 10 10 10 10 10 100 100 1000 1100 1200];
st2 = [ 10 10 100 101 101 1010 1090 1150 2000 3000];

matches = mysort.spiketrain.findNearSpikes(st1, st2, 0)

matches = mysort.spiketrain.findNearSpikes(st1, st2, 1)
% 
% matches = mysort.spiketrain.findNearSpikes(st1, st2, 10)
% 
% matches = mysort.spiketrain.findNearSpikes(st1, st2, 100)
% 
% matches = mysort.spiketrain.findNearSpikes(st1, st2, 1000)

st1 = sort(round(rand(1, 300000)*1000000));
st2 = st1(1:10:40000);

matches = mysort.spiketrain.findNearSpikes(st1, st2, 30);
