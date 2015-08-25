classdef preprocess < ana.botmpaper.harrisAnalysis.analysis
    properties
    end
    
    methods
        %------------------------------------------------------------------
        function self = preprocess()
            self = self@ana.botmpaper.harrisAnalysis.analysis('Preprocess');
        end
        %------------------------------------------------------------------
        function [fh figureNames] = makeFigures_(self, name, P)
            fh = []; p = 0; figureNames = {};
            R.I = self.getVariable(name, 'I');
            R.Xfil = self.getVariable(name, 'Xfil');
            R.Ifil = self.getVariable(name, 'Ifil');
            R.tgdf = self.getVariable(name, 'tgdf');
            R.gtSpikeTimes = self.getVariable(name, 'gtSpikeTimes');
            
            fh(p+1) = mysort.plot.figure('w', 1200, 'h', 1000); p=p+1;
                ah = subplot(2,1,1);
                plot(R.I(1:2:end));
                ah(2) = subplot(2,1,2);
                plot(R.Ifil(1:2:end));
                hold on
                plot(R.gtSpikeTimes/2, R.Ifil(R.gtSpikeTimes), 'rx', 'markersize', 20, 'linewidth', 3);
                linkaxes(ah, 'x');
                mysort.plot.figureTitle(name);
            figureNames{p} = 'IntracellularDetection';
            
            
            fh(p+1) = mysort.plot.figure('w', 1200, 'h', 1000); p=p+1;
                trueintraspikes = mysort.epoch.extractWaveform(R.I, [R.tgdf(:,2)-30 R.tgdf(:,2)+150]);
                trueextraspikes = mysort.epoch.extractWaveform(R.Xfil, [R.tgdf(:,2)-10 R.tgdf(:,2)+55]);
                subplot(2,1,1)
                plot(trueintraspikes')    
                subplot(2,1,2)
                plot(trueextraspikes')  
                mysort.plot.figureTitle(name);
            figureNames{p} = 'DetectedCutSpikes';
                
        end
    end
    methods(Access=protected)
        %------------------------------------------------------------------
        function result = runTrial_(self, name)
            [D X I] = ana.harris.util.load(name);
            
            hpf = 300; lpf= min(8000, D.srate/2); forder = 4; srate = D.srate;
            hd = mysort.mea.filter_design(hpf, lpf, srate, forder);
            Xfil = filtfilthd(hd, X')';

            hpf = 10; lpf= min(150, D.srate/2); forder = 4; srate = D.srate;
            hdI = mysort.mea.filter_design(hpf, lpf, srate, forder);
            Ifil = filtfilthd(hdI, I')';
            
            ifo = ana.harris.info(name);
            thr = ifo{10};
            [pks gtSpikeTimes] = findpeaks(Ifil+.0001*randn(size(Ifil)), 'minpeakheight', thr, 'minpeakdistance', 10);
            tgdf = [ones(length(gtSpikeTimes),1) gtSpikeTimes(:)-20];
            
            result.D = D;
            result.X = X;
            result.I = I;
            result.Xfil = Xfil;
            result.Ifil = Ifil;
            result.gtSpikeTimes = gtSpikeTimes;
            result.tgdf = tgdf;
        end
    end
end