k1 = [0 1 5];
k2 = {@length, @(x) size(x,1)};
k3 = {'peter', 'pan', 'paul'};
[p1 p2 p3] = mysort.util.combinationIterator(k1, k2, k3);