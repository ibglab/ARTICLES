function [ bases, nBases ] = generateEqualBoxcarBases( totalTime, binWidth)
% Generates boxcar basis function, each one of width specified in binWidth,
% wich together span temporal totalTime

nBases = round(totalTime / binWidth);
binIdx = [1:binWidth];
bases = zeros(nBases*binWidth, nBases);
for i = 1:nBases
    indexes = binIdx + binWidth*(i-1);
    bases(indexes, i) = 1;   
end

end

