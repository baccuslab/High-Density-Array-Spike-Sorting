% function testHdf5ReadWrite()

    R.foo = 'bar';
    R.C = {'a' 'ab'};
    R.D = {randn(2,10) randn(2,11)};
    R.a(1).bla = [1:10];
    R.a(2).bla = [1:15];
    R.empty = [];
    R.bss = mysort.sorters.BSSMP();
    R.bss.D = randn(2,1000);
    fname = 'bla.hd5';
    mysort.util.hdf5recursiveSave(fname, R);
    B = mysort.util.hdf5recursiveLoad(fname);
    
    