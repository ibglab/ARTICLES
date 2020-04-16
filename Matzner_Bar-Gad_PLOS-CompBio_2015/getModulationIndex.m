function [m] = getModulationIndex(peakPower, r0, T , winLen, fs)
% [m] = getModulationIndex(peak, r0, T , winLen, fs)
% This function takes as input the parameters of a spike train, and its
% spectrum (estimated using welch's method), and returns the modualtion index of the spike train.
% Input parameters:
%   peakPower: the spectral peak's power.
%   r0: base firing rate. (spike/s)
%   T: total recording time (length of the spike train).
%   winLen: the window length used in the calculation of the spectrum.
%   fs: the sampling rate of the spike train, in Hz.
%
% Output parameters:  
%   m : The modulation index of the spike train

% First correct the peak for comparison with Welch's estimator
r0 = r0/fs;
peak = correctWelchPeak(peakPower, winLen, T, r0, fs);
if (peak<r0)
    m =0;
else
    m = 2*(sqrt(peak-r0));
    m = m/(r0*sqrt(T));
    m = abs(m);
end

end