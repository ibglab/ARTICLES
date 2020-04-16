function [modulation, originalModulation, isSignificant] = getModIndexWithRefPer(spikeTrain, refPeriod, fs, freq, winLen, numSTDs)
% This function includes a correction procedure that corrects the
% modulation index to accomodate for the refractory period.
% Input parameters:
%   spikeTrain: a vector of point process spike train in the form [0, 1, 0, 0, 1, 0,...], where 1 represents the occurence of a spike.
%   refPeriod: length of the absolute refractory period, in ms
%   fs: the sampling rate of the spike train, in Hz.
%   freq: the frequncy of the spectral peak.
%   winLen: the window length to use in the estimation of the spectrum
%  numSTDs - STDs for detecting significant modulation. Default is 2 STDs.
%
% Output parameters:  
%   modulation : The corrected modulation index of the spike train
%   originalModulation : The original modulation index of the spike train, without correction.
%  isSignificant - 1 if the MI is significant

if (nargin < 7)   
    numSTDs  = 2;    
end
 isSignificant = 0 ;
 
startBand = freq-2;
endBand = freq+2;
totalTime = length(spikeTrain);
firingRate = sum(spikeTrain)/(totalTime/1000);
[spectrum, freqRange,snr, peakPower, peakFreq] = powerSpectrum(spikeTrain, fs, startBand, endBand, winLen);
originalModulation=  getModulationIndex(peakPower, firingRate,totalTime  , winLen, fs);

modulation = originalModulation;
if (originalModulation==0)
    return;
end
if (originalModulation > 1)
    return;
end

% calculate the firing rate without refractory period
totalTimeWithoutRef = totalTime - refPeriod*sum(spikeTrain);
meanRateWithoutRef  = (sum(spikeTrain)/(totalTimeWithoutRef /fs))/fs;
        
time=[1/fs:1/fs:totalTime/fs];

ex=1;
N= 100;
while ex
        rateFunc = (modulation*meanRateWithoutRef) *(cos(2*pi*freq*time))+meanRateWithoutRef;
        amps = zeros(1,N);
        for j = 1:N
            st= generatePoissonTrain(totalTime, rateFunc, refPeriod);              
            firingRate = sum(st)/(totalTime/1000);
            [spectrum, freqRange,snr, peakPower, peakFreq] = powerSpectrum(st, fs, startBand, endBand, winLen);

            amps(j) =  getModulationIndex(peakPower, firingRate,totalTime  , winLen, fs);
        end
                reconstructedMods = mean(amps);
        if (reconstructedMods >=originalModulation)
            ex=0;
        else 
            modulation = modulation+0.01;
        end
end 


% Find significance level
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
miZero=  getModulationIndex(mean(peakPowers) + numSTDs*std(peakPowers), mean(firingRates),totalTime  , winLen, fs);
if ( originalModulation > miZero)
    isSignificant = 1;
end
