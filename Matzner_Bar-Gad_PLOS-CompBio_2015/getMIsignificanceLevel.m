function [significanceLevel] = getMIsignificanceLevel(spikeTrain, refPeriod, fs, peakFreq,  winLen, numSTDs, N)
    %Get spike train and return the significance level for detecting
    %significant modulation index, acoording to the specified numSTDs(default is 2)    % 
 % Input parameters:
%   spikeTrain: a vector of point process spike train in the form [0, 1, 0, 0, 1, 0,...], where 1 represents the occurence of a spike.
%   refPeriod: length of the absolute refractory period, in ms
%   fs: the sampling rate of the spike train, in Hz.
%   peakFreq: the frequncy of the spectral peak.
%   winLen: the window length to use in the estimation of the spectrum
%   numSTDs - STDs for detecting significant modulation. Default is 2 STDs.
%   N - number of simulated spike train to generate. Default is 1000.
%
% Output parameters:  
%   modulation : The corrected modulation index of the spike train
%   originalModulation : The original modulation index of the spike train, without correction.


if (nargin < 6)
    numSTDs = 2;
    N  = 1000;    
end
if (nargin < 7)   
    N  = 1000;    
end

totalTime = length(spikeTrain);

totalTimeWithoutRef = totalTime - refPeriod*sum(spikeTrain);
meanRateWithoutRef  = (sum(spikeTrain)/(totalTimeWithoutRef /fs))/fs;
rateFunc = meanRateWithoutRef.*ones(1, totalTime);          % zero modulation 

firingRates = zeros(1, N);
numOfSpikes = zeros(1, N);
peakPowers = zeros(1, N);

 for k = 1:N
            aTrain = generatePoissonTrain(totalTime, rateFunc, refPeriod);            
            [spectrum, freqRange,snr, peakPower, freq] =  powerSpectrum(aTrain, fs, peakFreq, peakFreq, winLen);
            peakPowers(k) = peakPower;
            numOfSpikes(k) = sum(aTrain);
            firingRates(k) = sum(aTrain)/(totalTime/1000);
 end
[correctedModulationIndex, modulationIndex] = getModIndexFromPeakPower(mean(peakPowers) + numSTDs*std(peakPowers), refPeriod, fs, totalTime, mean(firingRates), mean(numOfSpikes), peakFreq, winLen);        
significanceLevel = correctedModulationIndex;
