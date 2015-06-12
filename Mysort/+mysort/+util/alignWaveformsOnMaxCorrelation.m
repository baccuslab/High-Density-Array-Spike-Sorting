
function [ali tau] = alignWaveformsOnMaxCorrelation(XI,varargin)
warning('This function is depricated! Use mysort.wf.* instead!');
error('This function does not work anymore!');
    P.C = [];
    P.trunc = 0;
    P.debug = 0;
    P = mysort.util.parseInputs(P,'alignWaveformsOnMaxCorrelation',varargin);

    [Tf nC nT] = size(XI);
    if isempty(P.C)
        P.C = eye(Tf); %*nC);
    end
    
    F = mysort.util.inverseFilters(XI,P.C);
    zeroOutFilterArtefacts = true;
    xvsf = mysort.util.calculateXIvsF(XI,F,nC,zeroOutFilterArtefacts);
    [M I] = max(xvsf,[], 1);
    I = squeeze(I);
    I = I -Tf;
    
    % find minimal shift possibility
    min_shift_sum = sum(I,2);
    [minshift, min_shift_idx] = min(abs(min_shift_sum));
    tau = I(min_shift_idx,:);
    
    
    XM = mysort.wf.t2m(XI);
    XM = mysort.util.shiftMCRows(XM, -tau, nC, P.trunc);
    
    ali  = mysort.wf.m2t(XM, nC);
    
    if P.debug
        FM = mysort.util.t2m(F);
        FM = mysort.util.shiftMCRows(FM, -tau, nC, P.trunc);
        aliF = mysort.util.m2t(FM, nC);
        
        disp(I);
        RES = mysort.plot.XIvsF(XI,F,'XIvsF',xvsf,'title',0,'axistight',1);        
        hold on
        RES = mysort.plot.XIvsF(ali,aliF,'title',0,'axistight',1,'figure',0,...
            'axesHandles',RES.axesHandles,'color','g','holdOn',1);

        [M I] = max(RES.XIvsF,[], 1);
        I = squeeze(I);
        I = I -Tf;
        disp(I);
    end

