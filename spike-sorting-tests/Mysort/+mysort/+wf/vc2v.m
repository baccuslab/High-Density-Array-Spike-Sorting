function V = vc2v(VC, Tf)
    % converts the vector-channel representation VC of the data into the
    % vector representation. Rows of VC correspond to channels and single
    % waveforms are concatenates. But all waveforms in VC have their first
    % channel waveform in row one. Single channel waveforms are of length
    % Tf, thus VC has dimensionality nC x (Tf*nS) where nC is the number of
    % channels, Tf is the length of a single channel waveform and nS is the
    % number of single channel waveforms
    %
    % V:  v1c1 v1c2 v1c3 ... v1c_nC
    %     ...
    %     vMc1 vMc2 vMc3 ... vMc_nC
    %
    % VC: v1c1   v2c1   v3c1   ... vMc_1
    %     ...
    %     v1c_nC v2c_nC v3c_nC ... vMc_nC
    %
    nC = size(VC,1);
    nS = size(VC,2)/Tf;
    assert(round(nS) == nS, 'number of columns of VC must be a multiple of Tf!');
    
    V = reshape(VC', Tf, nS, nC);
    V = permute(V, [1 3 2]);
    V = reshape(V, Tf*nC, nS)';