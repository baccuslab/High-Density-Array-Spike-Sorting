function defs = definitions()
    defs.testfile_location = '/net/bs-filesvr01/export/group/hierlemann/temporary/ijones/Spike_Sorting_Files/';
    defs.testfile_h5  = 'Trace_id843_2011-12-06T11_27_53_5.stream.h5';
    defs.testfile_ntk = 'Trace_id843_2011-12-06T11_27_53_5.stream.ntk';
    if ~isempty(strfind(computer, 'PCWIN'))
        defs.testfile_location = 'C:\LocalData\Ian\SpikeSortingTest\';
        defs.testfile_h5  = 'Trace_id843_2011-12-06T11_27_53_5.stream.h5';
        defs.testfile_ntk = '';
    end