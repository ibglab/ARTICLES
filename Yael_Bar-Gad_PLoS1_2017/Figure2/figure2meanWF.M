function []=figure2meanWF ()
rawFileName='F150818-0001E006.mat';
rawChannel='CRAW_006';
txtFileSortedLinear='F150818-0001E006.matCRAW_006FiltFilt.txt';

load (rawFileName);
rawSig=eval(rawChannel);

Fs=44000; 

% Linear (zero) filtering using filtfilt
HP =300; LP =  6000;
[bbp, abp] = butter(4, [HP LP]/(Fs/2), 'bandpass');
linearFilteredSig=filtfilt(bbp, abp, double(rawSig));

%WAVEFORM
x=csvread(txtFileSortedLinear,1);
linearWFtimeStamps=x(find(x(:,2)==1),3);
linearWFtimeStamps=round((linearWFtimeStamps)*Fs);
timeWinBefore=0.5*Fs/1000;
timeWinAfter=1.5*Fs/1000;
linearWF=zeros(length(linearWFtimeStamps),timeWinBefore+timeWinAfter+1);

for i=1:length (linearWFtimeStamps)
    linearWF(i,:)=linearFilteredSig((linearWFtimeStamps(i)-timeWinBefore):(linearWFtimeStamps(i)+timeWinAfter));
end

meanLinearWF=(mean(linearWF)*38.147/200)-20;

figure2WFreconstruct (meanLinearWF);
