function alignmentTest()
    L = 2000;
    St1 = {[rand(1,100)*L] [rand(1,100)*L] [rand(1,100)*L]};
    St2 = {St1{1}(1:2:end)+2 St1{2}(2:2:end)-3 [rand(1,100)*L]};
    R = mysort.spiketrain.align(St1, St2, 5, 5,5);
    mysort.plot.alignment(R, 'restrict2Time', [800 1200])
    mysort.plot.printEvaluationTable(R)