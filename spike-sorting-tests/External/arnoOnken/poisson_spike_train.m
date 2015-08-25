function train = poisson_spike_train(rate, duration)

if (nargin < 2)
    error('poisson_spike_train: usage poisson_spike_train(rate, duration)')
end

if (~isscalar(rate) || rate < 0)
    error('poisson_spike_train: rate must be a non-negative scalar')
end
if (~isscalar(duration) || duration < 0)
    error('poisson_spike_train: duration must be a non-negative scalar')
end

if duration == 0
    train = [];
else
    % Generate spike train <train_poiss> with poisson statistics
    train = [];
    t_next = -log(rand(1)) / rate;
    while t_next <= duration
        train = [train, t_next];
        t_next = t_next - log(rand(1)) / rate;
    end
end
