function [classes centers nClasses obj R] = gmmClustering(X, varargin)
    P.kmin = 1;
    P.kmax = 6;
    P.repeats = 40;
    P.method = 'AIC';
    P.prewhitened = false;
    P.homoscedastic = true;
    P.debugFun = [];
    P.fixedCov = [];
    P = mysort.util.parseInputs(P, 'gmmClusteringBIC', varargin);
    
    covtype = 'full';
    if P.prewhitened
        covtype = 'diagonal';
    end
    
%     if ~isempty(P.fixedCov)
%         covtype = 'full';
%     end
    
    R.options = statset('Display','off');
    R.fitting_objs = {};
    R.k = P.kmin:P.kmax;
    R.AIC = zeros(1, length(R.k));
    R.BIC = R.AIC;
    for i = 1:length(R.k)
        [R.fitting_objs{i} R.AIC(i) R.BIC(i)] = mysort.clustering.fit(X,R.k(i),...
            'Options', R.options, 'Replicates', P.repeats,...
            'SharedCov', P.homoscedastic, 'CovType', covtype, ...
            'FixedSigma', P.fixedCov);
%         R.fitting_objs{i} = gmdistribution.fit(X,R.k(i),...
%             'Options', R.options, 'Replicates', P.repeats,...
%             'SharedCov', P.homoscedastic, 'CovType', covtype);        
%         R.AIC(i) = R.fitting_objs{i}.AIC;
%         R.BIC(i) = R.fitting_objs{i}.BIC;
    end
    
    [R.minAIC, R.minAICidx] = min(R.AIC);
    R.nClusterAIC = R.k(R.minAICidx);
    
    [R.minBIC, R.minBICidx] = min(R.BIC);
    R.nClusterBIC = R.k(R.minBICidx);

    if strcmp(P.method, 'BIC')
        obj = R.fitting_objs{R.minBICidx};
    else
        obj = R.fitting_objs{R.minAICidx};
    end
    classes = cluster(obj, X);
    centers = obj.mu;
    
    if nargout>2
        cid = unique(classes);
        nClasses = zeros(1,length(cid));
        for i=1:length(cid)
            nClasses(i) = sum(classes==cid(i));
        end
    end