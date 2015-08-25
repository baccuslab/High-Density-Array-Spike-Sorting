function vectorColorTest()
    fh = mysort.plot.figure('w', 500, 'h', 800); hold on
    set(fh, 'Name', 'vectorColorTest');
    [c C] = mysort.plot.vectorColor(0);
    numofcol=size(C,1);

    for i=numofcol:-1:0
        plot([0 1], [i i], '-', 'color',mysort.plot.vectorColor(i),'lineWidth',3) 
        strs{numofcol-i+1} = ['i: ' num2str(i)];
    end

    legend(strs)