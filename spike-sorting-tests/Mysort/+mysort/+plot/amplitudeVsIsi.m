function amplitudeVsIsi(A, st, varargin)
    P.axesHandles = [];
    P.IDs = [];
    P.maxIsi = 10000;
    P.srate = [];
    P.d3 = 0;
    P.dIsi = 1;
    P.plotParas = {'kx'};
    P = mysort.util.parseInputs(P, varargin, 'error');
    if isempty(A) || isempty(st)
        return
    end
    if ~iscell(A)
        A = {A};
    end
    if ~iscell(st)
        st = {st};
    end
    nId = length(A);
    assert(nId == length(st), 'number of spike trains and number of amplitude trains must be identical!');
    
    if isempty(P.axesHandles)
        P.axesHandles = gca;
    end
    if ~isempty(P.srate)
        zlabels = '2nd isi [ms]';
        xlabels = 'isi [ms]';
    else
        zlabels = '2nd isi [sample]';
        xlabels = 'isi [sample]';        
    end
    if isempty(P.srate); P.srate = 1; end
    fac = P.srate/1000;
    for i = 1:nId
        t = st{i};
        if size(t,2) > size(t,1); t=t'; end
        a = A{i};
        if size(a,2) > size(a,1); a=a'; end
        AT = sortrows([t a]);
        if size(AT,1) < 2
            continue
        end
        isi = diff(AT(:,1));
        maxisi = max(isi);
        if ~P.d3 
            if P.dIsi == 1
                plot(P.axesHandles(i), isi/fac, AT(2:end, 2), P.plotParas{:});
                axis(P.axesHandles(i), 'tight');
            else
                plot(P.axesHandles(i), isi(1:end-1)/fac, AT(3:end, 2), P.plotParas{:});
                axis(P.axesHandles(i), 'tight');
            end
        else
            %isi2 = isi(2:end)+isi(1:end-1);
            plot3(P.axesHandles(i), isi(2:end)/fac, AT(3:end, 2), isi(1:end-1)/fac, P.plotParas{:});
            zlabel(P.axesHandles(i), zlabels)
            axis(P.axesHandles(i), 'tight');
            set(P.axesHandles(i), 'zlim', [0 min(maxisi, P.maxIsi)/fac]);
        end
        xlabel(P.axesHandles(i), xlabels)
        ylabel(P.axesHandles(i), 'ampl');
        set(P.axesHandles(i), 'xlim', [0 min(maxisi, P.maxIsi)/fac]);
    end
    
    