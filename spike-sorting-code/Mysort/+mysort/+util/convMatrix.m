
function Mf = convMatrix(M, f)
    % convolutes ever row of M with f
    Mf = zeros(size(M));
    for i=1:size(M,1)
        Mf(i,:) = mysort.util.conv(M(i,:), f, 'same');
    end