%%
if 1
    clear mea
    clear mysort.mea.CMOSMEA
end

fstr = 'Trace_id843_2011-12-06T11_27_53_5.stream.h5';
mea = mysort.mea.CMOSMEA(fstr);

[nT nC] = size(mea);


bigChunkPositions = [1 1000 10000 100000 1000000 10000000];
bigChunkLen = 1000;

smallChunkPositions = [1 10 100 400 900];
smallChunkLen       = 100;

epochs = [1 1
          1 10
          5 90
          7 1
          55 10
          500 30
          990 10];
          
nP = length(bigChunkPositions);
% figure;
% for bc = 1:nP
%     ah = subplot(nP, 1, bc);
%     X = mea(1:3, bigChunkPositions(bc):bigChunkPositions(bc)+bigChunkLen-1);
%     spacer = mysort.plot.mc(X', 'figure', 0, 'axesHandle', ah, 'color', {'k'});
%     hold on
%     for sc = 1:length(smallChunkPositions)
%         s1 = bigChunkPositions(bc)+smallChunkPositions(sc)-1;
%         X = mea(1:3, s1:s1+smallChunkLen-1);
%         mysort.plot.mc(X', 'figure', 0, 'axesHandle', ah, 'color', {mysort.plot.vectorColor(sc)}, ...
%             'sampleOffset', smallChunkPositions(sc)-1, 'spacer', spacer);
%     end
%     for e = 1:size(epochs,1)
%         s1 = bigChunkPositions(bc)+epochs(e,1)-1;
%         X = mea(1:3, s1:s1+epochs(e,2)-1);
%         mysort.plot.mc(X', 'figure', 0, 'axesHandle', ah,...
%             'color', {mysort.plot.vectorColor(e+length(smallChunkPositions))}, ...
%             'linewidth', 2,...
%             'lineStyle', ':',...
%             'sampleOffset', epochs(e,1)-1, 'spacer', spacer);
%     end
% end

figure;
for bc = 1:nP
    ah = subplot(nP, 1, bc);
    X = mea(bigChunkPositions(bc):bigChunkPositions(bc)+bigChunkLen-1, end-2:end);
    spacer = mysort.plot.mc(X', 'figure', 0, 'axesHandle', ah, 'color', {'k'});
    hold on
    for sc = 1:length(smallChunkPositions)
        s1 = bigChunkPositions(bc)+smallChunkPositions(sc)-1;
        X = mea(s1:s1+smallChunkLen-1, end-2:end);
        mysort.plot.mc(X', 'figure', 0, 'axesHandle', ah, 'color', {mysort.plot.vectorColor(sc)}, ...
            'sampleOffset', smallChunkPositions(sc)-1, 'spacer', spacer);
    end
    for e = 1:size(epochs,1)
        s1 = bigChunkPositions(bc)+epochs(e,1)-1;
        X = mea(s1:s1+epochs(e,2)-1, end-2:end);
        mysort.plot.mc(X', 'figure', 0, 'axesHandle', ah,...
            'color', {mysort.plot.vectorColor(e+length(smallChunkPositions))}, ...
            'linewidth', 2,...
            'lineStyle', ':',...
            'sampleOffset', epochs(e,1)-1, 'spacer', spacer);
    end
end

%%
spikeTimes = [1 1000 10000 100000];
cutLeft = 30;
cutRight = 30;
spikeTimes = spikeTimes+cutLeft;
spikeLen = cutLeft+cutRight+1;

figure;
ah = axes();
hold on
X1 = mea.getCutWaveforms(spikeTimes, cutLeft, cutRight);
X2 = mea.getCutWaveforms(spikeTimes(2:end), cutLeft, cutRight);

for i=1:size(X1,1)
    mysort.plot.mc(mysort.util.v2m(X1(i,:), nC), 'figure', 0, 'axesHandle', ah, 'color', {mysort.plot.vectorColor(i)}, ...
        'sampleOffset', i*(spikeLen+10), 'spacer', spacer, 'linewidth', 2, 'lineStyle', ':');
end
for i=1:size(X2,1)
    mysort.plot.mc(mysort.util.v2m(X2(i,:), nC), 'figure', 0, 'axesHandle', ah, 'color', {mysort.plot.vectorColor(i)}, ...
        'sampleOffset', (i+1)*(spikeLen+10), 'spacer', spacer, 'linewidth', 2, 'lineStyle', '--');
end

for i=1:length(spikeTimes)
    X = mea(spikeTimes(i)-cutLeft:spikeTimes(i)+cutRight, 1:92);
    mysort.plot.mc(X', 'figure', 0, 'axesHandle', ah, 'color', {mysort.plot.vectorColor(i)}, ...
            'sampleOffset', i*(spikeLen+10), 'spacer', spacer);
        
end    