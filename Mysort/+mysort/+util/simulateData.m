function X = simulateData(T, gdf, noise)
    % Simulates a piece of data of length noise containing the templates
    % (rows in T, multiple channels concatenated) at positions defined in 
    % gdf. gdf is a matrix with two columns and one row per spike. the
    % first row is the index into the row in T (thus the template) the
    % second is the startsample of the spike in X.
    
    nC = size(noise,1);
    L  = size(noise,2);
    
    if ndims(T) == 2
        T = mysort.wf.v2t(T, nC);
    end
    
    nT = size(T,3);
    Tf = size(T,1);
    nS = size(gdf,1);
    minId = min(gdf(:,1));
    assert(minId>0, 'the gdf has an invalid id! (<1)');
    maxId = max(gdf(:,1));
    assert(maxId<=nS, 'the gdf has an invalid id! (>nS)');
    
    maxTime = max(gdf(:,2));
    assert(maxTime <= L-Tf, 'A spike wont fit inside the noise!');
    
    X = zeros(nC, L);
    
    if size(gdf,2) == 2
        gdf = [gdf ones(size(gdf,1),1)];
    end
    
    for i=1:nS
        id = gdf(i,1);
        t  = gdf(i,2);
        a  = gdf(i,3);
        X = mysort.util.shiftSubtract(X, -a*squeeze(T(:,:,id))', t);
    end
    X = X + noise;