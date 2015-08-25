
function a = getAnnotationLayer(h)
%First turn show hidden handles on:
set(0,'showhiddenhandles','on')

% then look for all axes in the figure of choice:
h_all_axes = findall(h,'type','axes');

%then get the 'annotation layer' axes handle:
a = double(find(handle(h_all_axes),'-class','graph2d.annotationlayer'));

% and turn off the hidden handles:
set(0,'showhiddenhandles','off');