function [Ts tau] = tComputeSubsampleShiftedVersions(tT, upsample)
    % Computes the upsampled version of tT and shifts them (up)sample wise.
    % Then, each possible subsample shifted version is downsampled again
    % and returned in T. 
    % Upsampling is done with Matlabs resample function, but also
    % mysort.wf.mSincfun would be possible
    %
    % Input:
    %    tT       - 
    %    nC       -
    %    upsample - 
    %
    % Output:
    %    Ts - 4 dimensional tensor, time x channels x templates x subshifts
    %  tau  - 
    [Tf nC nT] = size(tT);
    downsample = 1;
    doNotCutTail = 1;
    tTup = mysort.wf.tResample(tT, upsample, downsample, doNotCutTail);
    
    Ts = zeros(Tf, nC, nT, upsample);
    Ts(:,:,:,1) = tT;
    tau = [0:upsample-1]/upsample;
    for t=1:nT
        % the first one is the original register shift, already copied
        for subtaui=2:upsample
            Ts(:,:,t,subtaui) = tTup(subtaui:upsample:end,:,t);
        end
    end
    