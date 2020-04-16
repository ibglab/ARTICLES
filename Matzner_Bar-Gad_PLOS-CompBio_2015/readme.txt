This code is used in the article: "Quantifing spike train oscillation: biases. distortions and solutions", published in PLOS computational biology, 2015.

V1.0 - Ayala Matzner 14/4/2015

The package file contains 6 Matlab functions and 1 test script which runs them:

testScript.m
Test code for simulating an inhomogeneous Poisson spike train using an oscillatory rate function, and extracting the modulation index of this spike train.

1. getModulationIndex.m
This function takes as input the parameters of a spike train, and its spectrum (estimated using welch's method), and returns the modualtion index of the spike train.
Input parameters:
  peakPower: the spectral peak's power.
  r0: base firing rate (spikes/s).
  T: total recording time (length of the spike train).
  winLen: the window length used in the calculation of the spectrum.
  fs: the sampling rate of the spike train, in Hz.

Output parameters:  
  m : The modulation index of the spike train

2. correctedPeak.m
This function scales the power estimated by Welch's method to the analytically calculated power
Input parameters:
  peak: the spectral peak's power.
  winLen: the window length used in the calculation of the spectrum.
  T: total recording time (length of the spike train).
  r0: base firing rate.
  fs: the sampling rate of the spike train, in Hz.

Output parameters:  
  correctedPeak : the corrected peak's power.

3. getModIndexWithRefPer.m
This function includes a correction procedure that corrects the modulation index to accomodate for the refractory period.
Input parameters:
  spikeTrain: point process spike train
  refPeriod: length of the absolute refractory period, in ms
  fs: the sampling rate of the spike train, in Hz.
  freq: the frequncy of the spectral peak.
  winLen: the window length to use in the estimation of the spectrum
  numSTDs - STDs for detecting significant modulation. Default is 2 STDs.

Output parameters:  
  modulation : The corrected modulation index of the spike train
  originalModulation : The original modulation index of the spike train, without correction.
  isSignificant - 1 if the MI is significant

4. powerSpectrum.m
This function estimates the power spectrum of a spike train using Welch's method. It return  the peak in a specific searched band.
 Input parameters:
  spikeTrain: point process spike train
  fs: the sampling rate of the spike train, in Hz.
  stratBand: the low limit of the searched band, in Hz.
  endBand: the high limit of the searched band, in Hz.
  winLen: the window length to use in pwelch

Output parameters:  
  spectrum : The spectrum of the spike train
  freqRange: corresponding vector of frequncies, in Hz.
  snr: the SNR of the peak relative to the 100-500Hz baseline.
  peakPower: the power of the highest peak in the searched band
  peakFreq: the frequency of the highest peak in the searched band

5. generatePoissonTrain.m
This function generates a poisson spike train using a rate function

6.extractRequiredRecordingDuration.m 
This function calculates the recording duration required to detect a spectral peak.
Input parameters:
  r0: base firing rate.
  m: the modulation index of the spike train.
  snr: the desired SNR
  winLen: the window length used in the calculation of the spectrum.

Output parameters:  
  requiredTime: the required recording duration


