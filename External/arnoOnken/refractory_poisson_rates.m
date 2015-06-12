function rates = refractory_poisson_rates(lambda, A, duration, recovery_rate, trials)

[N, nlambda] = size(A);

rates = zeros(N, trials);
trains = cell(nlambda, 1);

for itrials = 1:trials
    % Generate trains and count spikes
    for itrain = 1:nlambda
        trains{itrain} = poisson_spike_train(lambda(itrain), duration);
    end
    for iN = 1:N
        train = [];
        for itrain = 1:nlambda
            if A(iN, itrain) == 1
                train = [train, trains{itrain}];
            end
        end
        % Thinning for refractory period
        train = sort(train);
        n = length(train);
        if n == 0
            train_ref = [];
        else
            train_ref = train(1);
        end
        for i = 2:n
            if (1 - exp((train(i-1) - train(i)) / recovery_rate)) >= rand(1)
                train_ref = [train_ref, train(i)];
            end
        end
        rates(iN, itrials) = length(train_ref);
    end
end
