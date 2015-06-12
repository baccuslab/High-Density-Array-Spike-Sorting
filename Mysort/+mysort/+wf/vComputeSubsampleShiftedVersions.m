function [Ts tau] = vComputeSubsampleShiftedVersions(vT, nC, upsample)
    % Computes the upsampled version of T shifts them (up)sample wise.
    % Then, each possible subsample shifted version is downsampled again
    % and returned in T. 
    % Upsampling is done with Matlabs resample function, but also
    % mysort.wf.mSincfun would be possible
    %
    % Input:
    %    T        - 
    %    nC       -
    %    upsample - 
    %
    % Output:
    %    Ts - 4 dimensional tensor, time x channels x templates x subshifts
    %  tau  - 
    tT = mysort.wf.v2t(vT, nC);
    [tT tau] = mysort.wf.tComputeSubsampleShiftedVersions(tT, upsample);
    Ts = mysort.wf.t2v(tT);
    
    