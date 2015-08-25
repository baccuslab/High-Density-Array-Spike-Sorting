
function [T, N_per_template] = calculateClassMeans(X, c, use_median)
    classes = unique(c);
    N_per_template = zeros(1, length(classes));
    T = zeros(length(classes), size(X,2));
    if nargin == 3 && ~isempty(use_median) && any(use_median>0)
        for i=1:length(classes)
            T(i,:) = median(X(c==classes(i),:), 1);
            N_per_template(i) = sum(c==classes(i));
        end
    else
        for i=1:length(classes)
            T(i,:) = mean(X(c==classes(i),:), 1);
            N_per_template(i) = sum(c==classes(i));
        end
    end