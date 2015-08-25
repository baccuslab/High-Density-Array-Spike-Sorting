function [ali tau xvsf] = tAlignOnCorrelation(T,varargin)
    P.C = [];
    P.trunc = 0;
    P.debug = 0;
    P.absMax = false;
    P = mysort.util.parseInputs(P,'alignWaveformsOnMaxCorrelation',varargin);

    [Tf nC nT] = size(T);
    if isempty(P.C)
        P.C = eye(Tf*nC);
    end
    vT = mysort.wf.t2v(T);
    vF = vT/P.C;
%     F = mysort.wf.v2t(vF, nC);
    
    zeroOutFilterArtefacts = true;
    xvsf = mysort.util.calculateXIvsF(vT,vF,nC,zeroOutFilterArtefacts);
    if P.absMax
        [M I] = max(abs(xvsf),[], 1);
    else
        [M I] = max(xvsf,[], 1);
    end
    I = squeeze(I);
    I = I -Tf;
    
    % find minimal shift possibility
    min_shift_sum = sum(I,2);
    [minshift, min_shift_idx] = min(abs(min_shift_sum));
    tau = I(min_shift_idx,:);
    
    
    XM = mysort.wf.t2v(T);
    XM = mysort.util.shiftMCRows(XM, -tau, nC, P.trunc);
    
    ali  = mysort.wf.v2t(XM, nC);
    
    if P.debug
        F = mysort.wf.v2t(vF, nC);
        FM = mysort.wf.v2m(vF, nC);
        FM = mysort.util.shiftMCRows(FM, -tau, nC, P.trunc);
        aliF = mysort.util.m2t(FM, nC);
        
        disp(I);
        RES = mysort.plot.XIvsF(T,F,'TvsF',xvsf,'title',0,'axistight',1);        
        hold on
        RES = mysort.plot.XIvsF(ali,aliF,'title',0,'axistight',1,'figure',0,...
            'axesHandles',RES.axesHandles,'color','g','holdOn',1);

        [M I] = max(RES.XIvsF,[], 1);
        I = squeeze(I);
        I = I -Tf;
        disp(I);
    end

