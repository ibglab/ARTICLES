function requiredTime = extractRequiredRecordingDuration (r0, m, snr,winLen)
% This function calculates the recording duration required to detect a
% spectral peak.
% Input parameters:
%   r0: base firing rate.
%   m: the modulation index of the spike train.
%   snr: the desired SNR
%   winLen: the window length used in the calculation of the spectrum.
%
% Output parameters:  
%   requiredTime: the required recording duration

requiredTime = (16*snr^2) / (winLen*r0^2*m^4);