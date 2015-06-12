
function [sorting SO] = sort(X, samp_per_sec)
    % myosort.sort(X) performs spike sorting on the data in X.
    % INPUT:
    %  X   - either a matrix containing multichannel data (rows correspond
    %        single recording channels) or a dataObject
    %  smap_per_sec  (optional) - samples per second for the recording in X
    %
    % OUTPUT:
    %  sorting - a matrix with two columns. Every row corresponds to one 
    %            identified spike in the recordings in X. The first column
    %            indicates the neuron ID, the second the time point for the
    %            spike
    %  SO  - sorting object. A class containing detaild information about
    %        the spike sorting and functions to access those details and
    %        visualize them
    import mysort;
    if nargin == 1
        samp_per_sec = [];
    end
    X = dataObject(X, samp_per_sec);
    
    SO = sortingObject(X);
    
    sorting = SO.sort();