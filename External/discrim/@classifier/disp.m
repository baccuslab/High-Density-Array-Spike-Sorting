function c = disp(f)
%CLASSIFIER/DISP Display CLASSIFIER object.

%   Copyright (c) 1999 Michael Kiefte.

%   $Log$

n = f.counts;
g = length(n);
prior = f.prior;
if isempty(prior)
  prior = n./sum(n);
elseif length(prior) == 1
  prior = repmat(1/g, 1, g);
end
range = f.range;
p = size(range, 2);

fprintf('    prior: ')
if isempty(prior)
  fprintf('(observed distribution)\n');
elseif length(prior) == 1
  fprintf('(equal priors)\n');
else
  fprintf('[');
  if g < 10
    fprintf('%.4g ', prior(1:g-1));
    fprintf('%.4g]\n', prior(g));
  else
    fprintf('1x%d double]\n', g);
  end
end

fprintf('   counts: ');
if length(n) == 1
  fprintf('%d (total counts)\n', n);
elseif g < 17
  fprintf('[');
  fprintf('%d ', n(1:end-1));
  fprintf('%d]\n', n(end));
else
  fprintf('1x%d double]\n', g);
end

fprintf('    range: [2x%d double]\n', p);
fprintf('     nvar: %d\n', p);
fprintf('   nclass: %d\n', g);



