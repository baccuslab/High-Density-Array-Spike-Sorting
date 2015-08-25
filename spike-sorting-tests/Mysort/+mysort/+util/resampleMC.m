
function Y = resampleMC(X,up,down, dimension)
    if nargin < 4
        dimension = 2;
    end
    if ndims(X) == 2
        if dimension == 1
            nC = size(X,2);
            for i=1:nC
                Y(:,i) = resample(X(:,i),up,down);
                if i==1
                    Y = [Y zeros(size(Y,1), nC-1)];
                end
            end           
        else
            nC = size(X,1);
            for i=1:nC
                Y(i,:) = resample(X(i,:),up,down);
                if i==1
                    Y = [Y; zeros(nC-1, size(Y,2))];
                end
            end               
        end
    elseif ndims(X) == 3
        for i=1:size(X,2)
            for k=1:size(X,3)
                Y(:,i,k) = resample(squeeze(X(:,i,k)), up, down);
            end
        end
    else
        error('not implemented');
    end
       