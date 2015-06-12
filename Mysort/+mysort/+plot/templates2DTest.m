
nC = 4;
Tf = 3;
T = 40*[ 0 .5 0  0  1  0   0 .5 0  0 0 0
         0 1  0  0  .5 0   0  0 0   0 1 0
         .1 .2 .1  0  2  0   0 2.1 0  0 0 0
         0 .4 0  0  .9  0   0 .4 0  0 0 0];
T = T(:, :);
tT = mysort.wf.v2t(T, nC);
elPos = [ 0  0
         20 20
         20  0
          0 20];

mysort.plot.templates2D(tT, elPos);