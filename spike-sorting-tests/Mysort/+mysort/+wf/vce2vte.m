function vte = vce2vte(vce, nC)
    % converts the channel embedding of a vectorized multichannel waveform
    % to the time embedding representation.
    % Input:
    %   vce  -  matrix, each row is one multi channel waveform. The samples
    %           of this waveform are organized as
    %           vce = [T1C1 T1C2 ... T1CN T2C1 ... T2CN ... TKC1 ... TKC1]
    %           where T is the time index (samples) and C is the channel
    %           index
    %   nC   -  Number of channels of the waveform (N in the example above)
    % Output:
    %   vte  -  usual vectorized multi channel waveform representation:
    %           vte = [C1T1 C1T2 ... C1TK C2T1 ... C2TK ... CNT1 ... CNTK]
    % Time embedding means
    %
    % x_t = x(a1 a2 ... aN b1 b2 ...)
    %
    % Here an example for four channels (a to d) and 3 timelags
    %
    % time embedding (x):
    %      a1 a2 a3  b1 b2 b3  c1 c2 c3  d1 d2 d3
    % tau: |--|--|
    % nC:     |---------|---------|---------|
    %
    % channel embedding (y):
    %      a1 b1 c1 d1  a2 b2 c2 d2  a3 b3 c3 d3
    % nC:  |--|--|--|
    % tau:      |------------|------------|    
    if isempty(vce)
        vte = [];
        return
    end
    
    Tf = size(vce,2)/nC;
    assert(Tf==round(Tf), 'nC does not match vce !');
    nT = size(vce,1);
    
    vte = reshape(permute(reshape(vce', [nC Tf nT]), [2 1 3]), [nC*Tf nT])';