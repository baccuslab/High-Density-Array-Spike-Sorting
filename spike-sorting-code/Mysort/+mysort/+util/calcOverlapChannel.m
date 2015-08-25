
function dd = calcOverlapChannel(D, M, overlapTypes)
% Computes the BOTM discriminant functions for overlapping spikes if the
% single template discriminant functions are in D and the different
% overlapping types are in overlapTypes. overlapTypes is an array of
% structs each vontaining the identity of the two filters (f1 and f2) and
% the time difference between them inside the overlap (tau). Furthermore,
% it contains the relative shifts of the individual templates (tau_f1 for
% f1 and tau_f2 for f2). tau_f1 - tau_f2 must be tau! So f1 normally shifts
% to the right (positive values of tau), f2 to the left (negative values of
% tau). M is the confusion tensor containing the cross correlation
% functions of all filters to all templates.
    len = size(D,2);
    nF  = size(D,1); 
    nOverlapChannel = length(overlapTypes)*(nF*(nF-1)/2);
    Tf = (size(M,1)-1)/2 + 1;
    dd = nan(nOverlapChannel, len);
    
    for i=1:length(overlapTypes)
        O = overlapTypes(i);
        nothing = zeros(size(D(1,:)));
        D1  = mysort.util.shiftSubtract(nothing, -D(O.f1,:), O.tau_f1);
        D2  = mysort.util.shiftSubtract(nothing, -D(O.f2,:), O.tau_f2);
        D12 = M(Tf + O.tau, O.f2, O.f1);
        dd(i,:) =  D1 + D2 - D12;
    end  