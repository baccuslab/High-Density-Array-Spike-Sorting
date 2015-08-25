function SC = openSliderWithSorting(varargin)
    DS = {};
    spacer = [];
    count = 1;
    for i=1:2:nargin
        DS{count} = varargin{i};
        x = DS{count}(1:min(10000, size(DS{count},1)), :);
        spacer(count) = 4*std(x(:));
        SC = mysort.spiketrain.SpikeSortingContainer('bla', varargin{i+1}, 'wfDataSource', DS{count});
        DS{count}.addSpikeSorting(SC);
        count = count+1;
    end
    mysort.plot.SliderDataAxes(DS, 'channelSpacers', spacer);