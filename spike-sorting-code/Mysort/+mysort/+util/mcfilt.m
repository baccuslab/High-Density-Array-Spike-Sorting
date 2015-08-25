
function xx = mcfilt(X, F, varargin)
    % calculates the multichannel convolution between matrix
    % X (data) and F (filter). Rows represent channels
    assert(size(X,1) == size(F,1), 'X and F must have same number of channels!')
    if nargin == 2
        varargin = {'same'};
    end
    for i=1:size(X,1)
        %xx(i,:) = xcorr(X(i,:), X(i,:), varargin{:});
        %xx(i,:) = filter(fliplr(F(i,:)), 1, X(i,:), varargin{:});
        if i==1
            xx = mysort.util.conv(X(i,:), fliplr(F(i,:)), varargin{:});
        else
            % ignore channels where filter is all zero
            if any(F(i,:))
                xx = xx + mysort.util.conv(X(i,:), fliplr(F(i,:)), varargin{:});
            end
        end
    end
end 