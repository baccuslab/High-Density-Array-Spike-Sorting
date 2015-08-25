
function num = getNumberInStrangeBasis(i, basisCardinalities)
    basisCardinalities = basisCardinalities(:)';
    num = zeros(1, length(basisCardinalities));
    cp = cumprod(fliplr(basisCardinalities));
    wertigkeit = fliplr([1 cp]);
    wertigkeit(1) =  [];
    for k=1:length(basisCardinalities)
        num(k) = floor(i/wertigkeit(k));
        i = i-num(k)*wertigkeit(k);
    end
end