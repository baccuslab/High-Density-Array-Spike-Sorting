% The refractory period obeys the equation
% recovery_rate * (dr/dt) = rate - r
function train = refractory_poisson_spike_train(rate, duration, recovery_rate)

if (nargin < 3)
    error('refractory_poisson_spike_train: usage refractory_poisson_spike_train(rate, duration, recovery_rate)')
end

if (~isscalar(rate) || rate < 0)
    error('refractory_poisson_spike_train: rate must be a non-negative scalar')
end
if (~isscalar(duration) || duration < 0)
    error('refractory_poisson_spike_train: duration must be a non-negative scalar')
end
if (~isscalar(recovery_rate) || recovery_rate < 0)
    error('refractory_poisson_spike_train: recovery_rate must be a non-negative scalar')
end

if duration == 0
    train = [];
else
    % Generate spike train <train_poiss> with poisson statistics
    train_poiss = [];
    t_next = -log(rand(1)) / rate;
    while t_next <= duration
        train_poiss = [train_poiss, t_next];
        t_next = t_next - log(rand(1)) / rate;
    end

    if (recovery_rate > 0)
        % Thinning; results in train with refractory period
        n = length(train_poiss);
        if n == 0
            train = [];
        else
            train = train_poiss(1);
        end
        for i = 2:n
            if (1 - exp((train_poiss(i-1) - train_poiss(i)) / recovery_rate)) >= rand(1)
                train = [train, train_poiss(i)];
            end
        end
    end
end
