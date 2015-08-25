function getBurstIndexTest()
    st = [10 20 30 100 110 120 130 200 300 300 301 302 400];
    idx = mysort.util.getBurstIndex(st, 10);
    gtidx = [1 2 3 1 2 3 4 1 1 2 3 4 1];
    assertEqual(idx, gtidx);