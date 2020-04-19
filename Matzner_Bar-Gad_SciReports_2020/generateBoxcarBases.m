function [ bases, nBases ] = generateBoxcarBases( totalTime, x)
% Generates boxcar basis function, each one of width
% 1x, 2x, 3x, 4x... where x is 
% specified in binWidth,
% wich together span temporal totalTime

i = 1;
startIndex = 1;
while (startIndex <= totalTime)
    endIndex =  startIndex + x*i - 1;   
    bases(startIndex : endIndex , i ) = 1;
    i = i + 1;
    startIndex = endIndex + 1;
end
nBases = i-1;
bases = flipud(bases);
end

