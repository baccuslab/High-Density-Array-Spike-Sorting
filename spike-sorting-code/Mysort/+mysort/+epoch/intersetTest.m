function intersetTest()
e1 = [1 3
      8 15
      20 22
      100 200];
  
e2 = [2 2
      1 3
      7 9
      14 16
      19 23
      110 120
      130 140
      190 210];
e3 = mysort.epoch.intersect(e1,e2)

mysort.plot.figure();
mysort.plot.epochs(e1, 2, 'b', 'linewidth', 4);
hold on
mysort.plot.epochs(e2, 1, 'r', 'linewidth', 4);
mysort.plot.epochs(e3, 0, 'm', 'linewidth', 4);
set(gca, 'ylim', [-1 3]);
