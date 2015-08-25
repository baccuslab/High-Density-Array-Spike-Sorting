
function marker = markerForClass(i)
% Just returns colors
C = [  0  0  1      % 01 blue
       0  1  0      % 02 lime
       1  0  0      % 03 red
       .9  .9  0      % 04 yellow
       1  0  1      % 05 fuchsia
       0  1  1      % 06 aqua
      .5 .5  1      % 08
      .7 .3 .3      % 09
      .3 .7 .3      % 10
       0  0 .5      % 11 navy
       0 .5  0      % 12 green
      .5  0  0      % 13 maroon
       0 .5 .5      % 14 teal
      .5  0 .5      % 15 purple
      .5 .5  0      % 16 olive
      .2 .5 .8      % 17
      1  .1 .1
      .3 .3 .3
      0  0   0      % 20 black      
      .9 .2 .9
      .2 .9 .9
      .9 .9 .2
      .5 .5 .5      % 07 gray
      .4 .7 .3];
nColors = size(C,1);
markerIdx = floor(i/nColors)+1;

markers = {'.', 'x', 'o', 'd'};
marker = markers{markerIdx};