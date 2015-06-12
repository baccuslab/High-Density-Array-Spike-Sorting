
function [xcovs D autocovs autocov_norms NE smad] = noise_covariance_over_distance(X,x,y, varargin)
P.Tf = 100;
P.thr = 4.5;
P.maxLag = 10;
P.nNeighbors = 6;
P.maxDist = 40;
P = mysort.util.parseInputs(P, 'noise_covariance_over_distance', varargin);

nC = size(X,1);
NE = cell(1,nC);     % store noise epochs
smad = zeros(1, nC); % store standard deviations
%% 
Y = X;
for ch = 1:nC
    [s, ne] = ana.douglas.estimateSigma(X(ch,:), P.Tf, P.thr);
    smad(ch) = s;
    NE{ch} = ne;
    Y(ch,:) = Y(ch,:)/s;
end

%% Compute noise epochs common with n nearest electrodes
NEN = NE;
for ch = 1:nC
    % find k nearest electrodes
    neighbors = mysort.mea.nearestElectrodes(x, y, ch, P.nNeighbors);
    for n = 1:length(neighbors)
        NEN{ch} = mysort.epoch.intersect(NEN{ch}, NE{neighbors(n)});
    end
end

%%
nXc = (nC*(nC-1))/2;
D = zeros(1,nXc);count = 1;

autocovs = zeros(nC, 2*P.maxLag+1);
autocov_norms = zeros(nC,1);
xcovs = zeros(nXc, 2*P.maxLag+1);
IDX = zeros(nXc,2);

for ch1 = 1:nC
    p1 = [x(ch1) y(ch1)];
    NE{ch1} = mysort.epoch.removeShort(NE{ch1}, 400);
    [ac lc] = mysort.util.xcorr_in_epochs(Y(ch1,:), NE{ch1}, P.maxLag);
    autocovs(ch1,:) = ac;
    autocov_norms(ch1) = lc;
    
    for ch2 = ch1+1:nC 
        p2 = [x(ch2) y(ch2)];
        d = norm(p1-p2);
        IDX(count,:) = [ch1 ch2];
        D(count) = d; 
        if d < P.maxDist
            cNE = mysort.epoch.intersect(NE{ch1},NE{ch2});
%             cNE = cNE(1:3:end,:);
            [xc lc] = mysort.util.xcorr_in_epochs(Y([ch1 ch2],:), cNE, P.maxLag);
        else
            xc = zeros(1, 2*P.maxLag+1);
        end
        xcovs(count,:) = xc;   
        count=count+1;
        disp(count)
    end
end


