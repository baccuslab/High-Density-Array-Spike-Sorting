function [Y t] = downsample_for_plot(X, pixel)
   % resamples a given data piece to a smalle data piece that is suitable
   % for screen representation
   % input:
   %    x     - data piece
   %    pixel - number of pixels available for the plot
   % output:
   %    y     - downsampled data of length 2*pixel
   %    t     - timepoints in sample of the old x for every sample in y

   if exist('pixel', 'var') == false
       pixel = 1000;
   end

   L  = size(X,2);
   nC = size(X,1);
   if L <= 2*pixel
       % length of data is already small enough. No resampling
       y = x;
       t = 1:L;
       return
   end

   % init result
   Y = zeros(nC, 2*pixel);
   t = zeros(1, 2*pixel);
   % build window positions
   windows = round(linspace(1, L, pixel+1));
   %win_len_half = round((windows(2)-windows(1))/2);
   % keep only min and max value for every window position
   for i=1:length(windows)-1
       [mi imi] = min(X(:, windows(i):windows(i+1)-1), [], 2);
       [ma ima] = max(X(:, windows(i):windows(i+1)-1), [], 2);
       a = -1; b = 0;
       if imi>ima; a = 0; b=-1; end
       Y(:, 2*i+a) = mi;
       Y(:, 2*i+b) = ma;
       t(1, 2*i+a) = windows(i)+imi-1;
       t(1, 2*i+b) = windows(i)+ima-1; %+win_len_half;%
   end