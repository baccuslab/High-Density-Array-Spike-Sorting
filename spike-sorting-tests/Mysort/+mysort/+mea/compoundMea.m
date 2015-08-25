function compoundMea = compoundMea(fnamelist, varargin)
    meas = cell(1,length(fnamelist));
    for i=1:length(fnamelist)
        meas{i} = mysort.mea.CMOSMEA(fnamelist{i}, varargin{:});
    end
    tmpMea = mysort.ds.MultiSessionMatrix('MultiH4Mea', meas, meas{1}.getSamplesPerSecond());
    tmpMea.concatenateSessions();
    compoundMea = mysort.ds.Matrix(tmpMea, meas{1}.getSamplesPerSecond());
    compoundMea.MultiElectrode = meas{1}.MultiElectrode;