function heatmaps(data, varargin)
P.fh = [];
P.nBins = [50 50];
if nargin == 2
    P.nBins = varargin{1};
else
    P = mysort.util.parseInputs(P, varargin, 'error');
end

if length(P.nBins) == 1
    P.nBins(2) = P.nBins;
end

nP = floor(size(data,2)/2);
if isempty(P.fh)
    figure
end

for i=1:nP
    i1 = 2*(i-1)+1;
    i2 = 2*(i-1)+2;
    [nc edgesx edgesy] = mysort.util.hist2d(data(:,i1:i2), P.nBins(:,1), P.nBins(:,2));
    subplot(2,nP,i);
    plot(data(:,i1), data(:,i2), '.');
    axis tight
    
    subplot(2,nP,nP+i);
    nc(nc==0) = nan;
    mysort.plot.imagesc(nc(1:end-1,1:end-1), 'xdata', edgesx(1:end-1), 'ydata', edgesy(1:end-1))
    axis normal xy
    colorbar
end
   
