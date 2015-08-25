
function [colvec C marker] = vectorColor(i)
    markerSet = mysort.plot.lineTypes();
    % Just returns colors
    C = [  0  0  1      % 01 blue
           0  1  0      % 02 lime
           1  0  0      % 03 red
           .9  .9  0    % 04 yellow
           1  0  1      % 05 fuchsia
           0  1  1      % 06 aqua
          .5 .5  1      % 07
          1 .5 .5 %.6 .3 .4      % 08
          .3 .7 .3      % 09
           0  0 .5      % 10 navy
           0 .5  0      % 11 green
          .5  0  0      % 12 maroon
           0 .5 .5      % 13 teal
          .5  0 .5      % 14 purple
          .5 .5  0      % 15 olive
          .2 .5 .8      % 16
          .8  .3 .2
          .3 .3 .3
          0  0   0      % 19 black      
          .7 .3 .8
          .2 .9 .9
          .8 .8 .3
          .5 .5 .5      % 07 gray
          .1 .3 .1];
      cidx = mod(i, size(C,1))+1;
      colvec = C(cidx,:);
      midx = min(length(markerSet), floor(i/size(C,1))+1);
      marker = markerSet{midx};
      
