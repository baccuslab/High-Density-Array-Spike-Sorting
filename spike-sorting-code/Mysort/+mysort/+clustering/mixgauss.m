
function [mu, SIG, p, LL,Z, HC, HCP, HCD]=mixgauss(x,k,varargin)
[N dim] = size(x);

%Default Parameters
p = ones(1,k)/k;      % mixing proportions
muinit='orig';
siginit='var';
SIG = zeros(dim,dim,k);  % squared covariance matrices
maxiter=200;            % number of iterations
verbose=false;
tol=1e-4;


for i = 1:2:nargin-2                   
    switch varargin{i}
        case 'muinit', muinit= varargin{i+1};
        case 'siginit',siginit=varargin{i+1};
        case 'tol', tol= varargin{i+1};
        case 'graph3d',graph3d=varargin{i+1};
        case 'output',
            if strcmp(varargin{i+1},'verbose'), verbose=true;
            elseif strcmp(varargin{i+1},'brief'), verbose=false;
            end
    end
end
%graph3d=graph3d & q==3 & k>1;
if strcmp(muinit,'rand'),
    mu = randn(dim,k);  
elseif strcmp(muinit,'orig')
    mu = zeros(dim,k);
else
    mu = muinit';
end

% initialize
for i=1:k
    if strcmp(siginit,'rand'),
        r=0.3*abs(rand(q,1));
        SIG(:,:,i) =cov(x,1)+r*r';
        %SIG(:,:,i) = -100*diag(log(rand(q,1))); 
    else
        SIG = siginit;
    end

end

LL=0;
finished=false;
t=0;
if graph3d,
    figure;
    title('MOG');
end
while ~finished,
    t=t+1;
    %fprintf('t=%d ',t);
    %E-step: Evaluate the responsibilities using the current parameter
    %values
   
    for i=1:k,
        %Z(:,i) = p(i)*det(SIG(:,:,i))^(-0.5)*exp(-0.5*sum((x'-repmat(mu(:,i),1,N))'*inv(SIG(:,:,i)).*(x'-repmat(mu(:,i),1,N))',2));
        Z(:,i) = p(i)*det(SIG)^(-0.5)*exp(-0.5*sum((x'-repmat(mu(:,i),1,N))'*inv(SIG).*(x'-repmat(mu(:,i),1,N))',2));
    end
    normalizer = repmat(sum(Z,2),1,k);
    normalizer(normalizer==0) = 100*eps;
    Z = Z./normalizer;
    %[ma hardClust] = max(Z,[],2);
    % Do the M-step: Re-estimate the parameters using the current responsibilities
    L=0;
    for i=1:k,
%         if t>1
%             myIdx = hardClust == i;
%             if ceil(100*sum(myIdx)/N)<5
%                 mu(:,i) = x(ceil(N*rand(1)),:); 
%                 L = 0;
%                 break
%             end
%         end
        
        % calculate the sample covariance matrix of component i
        Nk=sum(Z(:,i));

        mu(:,i) = (x'*Z(:,i))./Nk;
        
%         DISTS = sum( (x-repmat(mu(:,i)',N,1)).^2 ,2); 
%         [mi I] = min(DISTS);
%         mu(:,i) = x(I,:);

        p(i) = mean(Z(:,i));
        %SIG(:,:,i)=(x'-repmat(mu(:,i),1,N))*(repmat(Z(:,i),1,q).*(x'-repmat(mu(:,i),1,N))')./sum(Z(:,i));

        %L=L+p(i)*(det(2*pi*SIG(:,:,i))^(-0.5)*exp(-0.5*sum((x'-repmat(mu(:,i),1,N))'*inv(SIG(:,:,i)).*(x'-repmat(mu(:,i),1,N))',2)));    
        L=L+p(i)*(det(2*pi*SIG)^(-0.5)*exp(-0.5*sum((x'-repmat(mu(:,i),1,N))'*inv(SIG).*(x'-repmat(mu(:,i),1,N))',2)));    
    end
    LL(t)=sum(log(L));
    %LL(t)=log(sum(exp(loglik)));
    %fprintf('LL: %f \n',LL(t));

    if t>1,
        if abs(LL(t)-LL(t-1))<tol || t>maxiter,
            finished=true;
        end
    end
end
mu = mu';
[HCP HC] = max(Z,[],2);
CC = zeros(k, dim);
HCD = 0;
for i=1:k
    myIDX = HC==i;
    myN = sum(myIDX);
    if myN>0
        CC(i,:) = mean(x(myIDX,:));
        HCD = HCD + sum(mysort.util.maha(x(myIDX,:), repmat( CC(i,:),myN,1 ), SIG));
    end
end
for i=1:k
    myIDX = HC==i;
    myN = sum(myIDX);
    for j=1:k
        if j~=i
            if max(abs(CC(j,:)))>0
                HCD = HCD + sum(mvnpdf(x(myIDX,:), repmat( CC(j,:),myN,1 ), SIG));
            end
        end
    end
end


% -------------------------------------------------------------------------