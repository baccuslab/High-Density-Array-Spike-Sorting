function save2HDF5File(fname, sessionindex, sig, channel_list, samplesPerSecond, gain)
    if isempty(gain)
        gain = [958.558 30.03 31.92 1.0]';
    end
    S.sig = sig;
    
    channel_list = int32(channel_list);
    names = {'channel_nr', 'connected', 'x', 'y', 'idx', 'dummy', 'damaged'};
    for i=1:size(channel_list,1)
        CL(i) = hdf5.h5compound();
        for n=1:length(names)
            CL(i).addMember(names{n});
            set(CL(i), 'Data', {channel_list(i,1), channel_list(i,2), channel_list(i,3), channel_list(i,4), channel_list(i,5), channel_list(i,6), channel_list(i,7) })
        end
    end    
    S.channel_list = CL;
    S.gain = gain;
    S.sr = samplesPerSecond;
    S.version = .1;
    S.message = 'no message';
    S.chipid = 0;
    
    h5path = ['/Sessions/Session' num2str(sessionindex)];
    
    mysort.h5.recursiveSave(fname, S, h5path, 'append');