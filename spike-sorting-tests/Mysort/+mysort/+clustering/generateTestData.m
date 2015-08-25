function [X P T] = generateTestData(varargin)
    P.k = [];
    P.C = [];
    %P.homoscedastic = [];
    P.N = [];
    P.dim = 3;
    P.mu = [];
    P.mu_spread = 40;
    %P.uniform_cluster_fraction = 0;
    P = mysort.util.parseInputs(P, '', varargin);
    
    %% Init
    % Check number of clusters
    if isempty(P.k)
        P.k = ceil(rand*10);
    end
    
    % Check number of elements for each cluster
    if isempty(P.N)
        P.N = zeros(1,P.k);
        for i=1:P.k
            P.N(i) = ceil(rand*300)+50;
        end
    end

    % Check covariance matrix
    if isempty(P.C)
        P.C = eye(P.dim);
    end
    
    % Check cluster centers
    if isempty(P.mu)
        C_mu_spread = P.mu_spread*eye(P.dim);
        P.mu = mvnrnd(zeros(1, P.dim), C_mu_spread, P.k);
    end
    
    % generate training data
    X = zeros(sum(P.N), P.dim);
    P.classes = zeros(sum(P.N), 1);
    for i=1:P.k
        idx = sum(P.N(1:i-1))+1:sum(P.N(1:i));
        X(idx,:) = ...
            mvnrnd(P.mu(i,:), P.C, P.N(i));
        P.classes(idx) = ones(P.N(i),1)*i;
    end    
    
    if nargout == 3
        % generate test data
        T = zeros(sum(P.N), P.dim);
        for i=1:P.k
            idx = sum(P.N(1:i-1))+1:sum(P.N(1:i));
            X(idx,:) = ...
                mvnrnd(P.mu(i,:), P.C, P.N(i));
        end    
    end