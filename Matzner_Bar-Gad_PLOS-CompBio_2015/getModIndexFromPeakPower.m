function [modulation, originalModulation] = getModIndexFromPeakPower(peakPower, refPeriod, fs, totalTime, firingRate, nSpikes, freq, winLen)
% Function used for getMIsignificanceLevel.m

modIndex = getModulationIndex(peakPower, firingRate, totalTime , winLen, fs);
originalModulation = modIndex;

modulation = modIndex;
if (modIndex==0)
    return;
end
if (modIndex > 1)
    return;
end

  totalTimeWithoutRef = totalTime - refPeriod*nSpikes;
  meanRateWithoutRef  = (nSpikes/(totalTimeWithoutRef /fs))/fs;    
        
    time=[0:1/fs:totalTime/fs];
    time = time(1:end-1);
    ex=1;
    N= 100;
    while ex
        rateFunc = (modulation*meanRateWithoutRef) *(cos(2*pi*freq*time))+meanRateWithoutRef;       
        peakPowers = zeros(1,N);
        firingRates = zeros(1,N);
        for j = 1:N
            st= generatePoissonTrain(totalTime, rateFunc, refPeriod);              
            firingRates(j) = sum(st)/(totalTime/1000);
            [spectrum, freqRange,peakStd, peakPowers(j), peakFreq] = powerSpectrum(st, fs, freq, freq, winLen); 
        end
        reconstructedMods =  getModulationIndex(mean(peakPowers), mean(firingRates),totalTime  , winLen, fs);
        if (reconstructedMods >=modIndex)
            ex=0;
        else 
            modulation = modulation+0.01;
        end
    end 

