function filter_init(f, X)
    reset(f);
	f.PersistentMemory=1;
%     m = mean( X, 1);   %first 50ms
%     f.States = m([1 1],:); %2nd order bandpass (FIXME for others...)
    y = filter(f, X);
    y = filter(f, X);