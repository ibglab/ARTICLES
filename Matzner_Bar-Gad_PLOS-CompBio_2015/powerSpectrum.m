function [spectrum, freqRange,snr, peakPower, peakFreq] = powerSpectrum(spikeTrain, fs, startBand, endBand, winLen)
% This function estimates the power spectrum of a spike train using Welch's
% method. It return  the peak in a specific searched band.
%  Input parameters:
%   spikeTrain: point process spike train
%   fs: the sampling rate of the spike train, in Hz.
%   stratBand: the low limit of the searched band, in Hz.
%   endBand: the high limit of the searched band, in Hz.
%   winLen: the window length to use in pwelch
%
% Output parameters:  
%   spectrum : The spectrum of the spike train
%   freqRange: corresponding vector of frequncies, in Hz.
%   snr: the SNR of the peak relative to the 100-500Hz baseline.
%   peakPower: the power of the highest peak in the searched band
%   peakFreq: the frequency of the highest peak in the searched band

overlap=0;     
nfft = winLen;
[spectrum, freqRange] = pwelch(spikeTrain,winLen,overlap,nfft,fs);

startBand = find(freqRange <=startBand, 1, 'last');
endBand = find(freqRange >=endBand, 1, 'first');

 % check for significance
    meanFreqStart = 100;
    meanFreqStart = find(freqRange>=meanFreqStart,1,'first');
    meanFreqEnd = 500;
    meanFreqEnd = find(freqRange>=meanFreqEnd,1,'first');
    meanPower = mean(spectrum(meanFreqStart:meanFreqEnd));
    stdPower = std(spectrum(meanFreqStart:meanFreqEnd));
    [peakPower,index] = max(spectrum(startBand:endBand));     % look for peak in specific band
    peakFreq = freqRange(startBand+index-1);
    snr = (peakPower-meanPower) / stdPower;
    
end