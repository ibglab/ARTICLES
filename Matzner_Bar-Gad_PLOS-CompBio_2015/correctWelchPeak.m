function correctedPeak = correctWelchPeak(peak, winLen, T, r0,fs)
% This function scales the power estimated by Welch's method to the
% analytically calculated power
% Input parameters:
%   peak: the spectral peak's power.
%   winLen: the window length used in the calculation of the spectrum.
%   T: total recording time (length of the spike train).
%   r0: base firing rate.
%   fs: the sampling rate of the spike train, in Hz.
%
% Output parameters:  
%   correctedPeak : the corrected peak's power.

w= hamming(winLen);
U = w'*w;
welchFactor = (mean(w)*winLen)/U;
nWins = ceil(T/winLen);

peak = peak*fs/2;
peak = (nWins*peak - ((nWins-1)*r0));
correctedPeak = peak*welchFactor; 
