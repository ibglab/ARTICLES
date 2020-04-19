function [] = fig1A()

rawFileName='F150818-0001E006.mat';
rawChannel='CRAW_006';
txtFileSortedLinear='F150818-0001E006.matCRAW_006FiltFilt.txt';

load (rawFileName);
rawSig=eval(rawChannel);

sampRate = 44000;                   % Samples per second
msRate = sampRate/1000;   % Samples per millisecond 
micVCorr = 200/38.147;          % Correction to microVolt units

figure; hold on;

% Plot the raw (wide band pass) signal
subplot (3,4,2:3);
plot ([0:length(rawSig)-1]/sampRate, rawSig/micVCorr, 'k');
xlabel ('Time (s)');
ylabel ('V (µV)');
title ('Raw Signal');

% Linear (zero) phase using Butterwoth filtfilt
HP =300; LP =  6000;
[bbp, abp] = butter(4, [HP LP]/(sampRate/2), 'bandpass');
linearFilteredSig=filtfilt(bbp, abp, double(rawSig));
axis ([1.55 1.65 -300 200]);

subplot (3,4,7:8);
plot ([0:length(linearFilteredSig)-1]/sampRate, linearFilteredSig/micVCorr, 'k');
xlabel ('Time (s)');
ylabel ('V (µV)');
title ('Linear Phase (LP) Filter ');
axis ([1.55 1.65 -250 250]);

% Non linear phase using Butterwoth filter
[b, a] = butter(4, [HP LP]/(sampRate/2), 'bandpass');
nonlinearFilteredSig=filter(b, a, double(rawSig));
subplot (3,4,5:6);
plot ([0:length(nonlinearFilteredSig)-1]/sampRate, nonlinearFilteredSig/micVCorr, 'k');
xlabel ('Time (s)');
ylabel ('V (µV)');
title ('Non Linear Phase (NLP) Filter ');
axis ([1.55 1.65 -250 250]);

%WAVEFORM
x=csvread(txtFileSortedLinear,1);
linearWFtimeStamps=x(find(x(:,2)==1),3);
linearWFtimeStamps=round((linearWFtimeStamps)*sampRate);
timeWinBefore = 1*msRate;
timeWinAfter = 1.5*msRate;
rawWF=zeros(length(linearWFtimeStamps),timeWinBefore+timeWinAfter+1);
linearWF=zeros(length(linearWFtimeStamps),timeWinBefore+timeWinAfter+1);
nonlinearWF=zeros(length(linearWFtimeStamps),timeWinBefore+timeWinAfter+1);

for i=1:length (linearWFtimeStamps)
    rawWF(i,:)=rawSig((linearWFtimeStamps(i)-timeWinBefore):(linearWFtimeStamps(i)+timeWinAfter));
    linearWF(i,:)=linearFilteredSig((linearWFtimeStamps(i)-timeWinBefore):(linearWFtimeStamps(i)+timeWinAfter));
    nonlinearWF(i,:)=nonlinearFilteredSig(round((linearWFtimeStamps(i)+0.05*44-timeWinBefore)):round((linearWFtimeStamps(i)+0.05*44+timeWinAfter)));
end
subplot (3,4,10:11); hold on;
meanRawWF=(mean(rawWF)/micVCorr)-10;
meanLinearWF=(mean(linearWF)/micVCorr)-20;
stdLinearWF=std(linearWF)/micVCorr;
meanNonLinearWF=mean(nonlinearWF)/micVCorr;
stdNonLinearWF=std(nonlinearWF)/micVCorr;
plot([0:size(meanRawWF,2)-1]/msRate,meanRawWF,'k');
xL = [0:size(meanLinearWF,2)-1]/msRate;
plot(xL,(meanLinearWF),'r');
plot(xL,(meanLinearWF+stdLinearWF),'r--');
plot(xL,(meanLinearWF-stdLinearWF),'r--');
xNL = [0:size(meanNonLinearWF,2)-1]/msRate;
plot(xNL,(meanNonLinearWF),'b');
plot(xNL,(meanNonLinearWF+stdNonLinearWF),'b--');
plot(xNL,(meanNonLinearWF-stdNonLinearWF),'b--');
title ('Waveform');
xlabel('Time (ms)');
ylabel ('V (µV)');
axis tight;



