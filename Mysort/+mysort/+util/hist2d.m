function [nc edgesx edgesy] = hist2d(data, edgesx, edgesy)

if nargin < 2 || isempty(edgesx)
    edgesx = 100;
end

if nargin < 3 || isempty(edgesy)
    edgesy = 100;
end

if length(edgesx) == 1
    edgesx = linspace(min(data(:,1)), max(data(:,1)), edgesx);
    edgesx(end) = edgesx(end) +1;
end

if length(edgesy) == 1
    edgesy = linspace(min(data(:,2)), max(data(:,2)), edgesy);
    edgesy(end) = edgesy(end) +1;
end
 

xIn = zeros(length(edgesx)-1, size(data,1));
for i=1:length(edgesx)-1
    xIn(i,:) = data(:,1) >= edgesx(i) & data(:,1) < edgesx(i+1);
end
yIn = zeros(length(edgesy)-1, size(data,1));
for i=1:length(edgesy)-1
    yIn(i,:) = data(:,2) >= edgesy(i) & data(:,2) < edgesy(i+1);
end

nc = zeros(length(edgesy), length(edgesx));
for i=1:length(edgesx)-1
    for j=1:length(edgesy)-1
        nc(j, i) = sum(xIn(i,:) & yIn(j,:));
    end
end