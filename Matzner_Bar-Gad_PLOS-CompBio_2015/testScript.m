% testscriptm
%
% Test code for simulating an inhomogeneous Poisson spike train using an
% oscillatory rate function, and extracting the modulation index of this
% spike train,

% simulate  1 second of a 12Hz oscillatory spike train.
totalTime = 60000;
freq = 12;
fs = 1000;
time=[0:1/fs:totalTime/fs];
baseRate = 0.03;
modulation = 0.25;
refPeriod = 0;
winLen = 1000;
rateFunc = (modulation*baseRate) * (cos(2*pi*freq*time))+baseRate;
spikeTrain = generatePoissonTrain(totalTime, rateFunc, refPeriod);
firingRate = sum(spikeTrain)/(totalTime/1000);
[spectrum, freqRange,snr, peakPower, peakFreq] = powerSpectrum(spikeTrain, fs, freq-2, freq+2, winLen);
modulationIndex = getModulationIndex(peakPower,firingRate, totalTime , winLen, fs);

disp('Rate function modulation:');
disp(modulation);
disp('Extracted modulation index:');
disp(modulationIndex);


% simulate spike train with refractory period:
refPeriod = 2;
spikeTrain = generatePoissonTrain(totalTime, rateFunc, refPeriod);
[spectrum, freqRange,snr, peakPower, peakFreq] = powerSpectrum(spikeTrain, fs, freq-2, freq+2, winLen);
[correctedModulationIndex, originalModulationIndex, isSignificant] = getModIndexWithRefPer(spikeTrain, refPeriod, fs, peakFreq, winLen);

disp('Extracted Modulation index with ref. priod:');
disp(originalModulationIndex);
disp('Corrected Modulation index with ref. priod:');
disp(correctedModulationIndex);

