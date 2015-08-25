function b = isEventMap(X)
    b = isstruct(X) && isfield(X, 'event_map');