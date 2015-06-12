% Plots the interspike interval histogram of the spike train <train> for
% bins of size <bin_size>. train can be a 2-dimensional cell array with
% trains als elements.
% e.g interspike_interval_histogram(refractory_poisson_spike_train(86, 500, 0.005), 0.001)
function interspike_interval_histogram(train, bin_size)

isi = [];
if iscell(train)
    % We have several trains
    s = size(train);
    for i = 1:s(1)
        for j = 1:s(2)
            ctrain = train{i, j};
            nt = length(ctrain);
            for k = 2:nt
                isi = [isi, ctrain(k) - ctrain(k - 1)];
            end
        end
    end
else
    nt = length(train);
    for i = 2:nt
        isi = [isi, train(i) - train(i - 1)];
    end
end

max_isi = max(isi);
x = 0:bin_size:max_isi;
n = histc(isi, x);
% Transform to percent
n = n ./ nt .* 100;
bar(x, n);
xlabel('Interspike interval [s]');
ylabel('Percent of ISI');
xlim([0 0.1]);