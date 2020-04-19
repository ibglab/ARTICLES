function []= figure4C ()

Fs=44000;
T=120; %s

%SIMULATED DATA
%waveform
load('D:\Dorin\filters paper\IC data\STR_AP_3.mat');
samplingRateIC=20000;
ICspk=c011_AxoPatch(3650:3800);
derivativeIC=-diff(ICspk);

spkWF=interp (derivativeIC, round(Fs*10/samplingRateIC));
spkWF=resample(spkWF,1,10);
spkWF=spkWF*8.5;

%spkTime
spkTrain=poissonGenerator (T*1000, 15, 1, 1);
spkTime=find(spkTrain==1)./1000; %s

%spk signal
spkSig=zeros(1,T*Fs);
x=round(spkTime*Fs);
for i=1:length (x)
   spkSig ((x(i)-80):(x(i)-80+length(spkWF)-1))=spkSig ((x(i)-80):(x(i)-80+length(spkWF)-1))+spkWF;    
end

% noise
pinkNoise=pinknoise (T*Fs);
pinkNoise=(pinkNoise./40);

signal=(pinkNoise+spkSig)*1000;


% NLP filter
HP =300;
LP =  6000;
[b, a] = butter(4, [HP LP]/(Fs/2), 'bandpass');
nonlinearFilteredSig=filter(b, a, signal);

%CORRECTED
correctedSignal= fliplr(filter(b,a, fliplr(nonlinearFilteredSig)));

%meanWF
for i=1:length (spkTime)
    rawWF(i,:)=signal((x(i)-80):(x(i)-80+length(spkWF)-1));
    nonlinearWF(i,:)=nonlinearFilteredSig((x(i)-80):(x(i)-80+length(spkWF)-1));
    correctedWF(i,:)=correctedSignal((x(i)-80):(x(i)-80+length(spkWF)-1));
end

meanRawWF=(mean(rawWF));
meanNonLinearWF=mean(nonlinearWF);
meanCorrectedWF=mean(correctedWF)-0.015;

figure; 
xmin=50;
xmax=50.2;
ymin=-400;
ymax=400;

subplot(2,3,[1:2]);hold on;
plot ([0:length(signal)-1]/Fs,signal,'k');
plot ([0:length(nonlinearFilteredSig)-1]/Fs,nonlinearFilteredSig,'b');
axis ([xmin xmax ymin ymax]);
plot ([0:length(correctedSignal)-1]/Fs,correctedSignal,'g');

xlabel('Time (s)');
ylabel ('V (µV)');
title ('Simulated data');


subplot(2,3,3);hold on;
plot([0:size(meanRawWF,2)-1]/(Fs/1000),meanRawWF,'k');
plot([0:size(meanNonLinearWF,2)-1]/(Fs/1000),(meanNonLinearWF),'b');
plot([0:size(meanCorrectedWF,2)-1]/(Fs/1000),(meanCorrectedWF),'g');
v=axis;
axis ([1 3.5 v(3) v(4)]);
xlabel('Time (ms)');
ylabel ('V (µV)');

clear;


%REAL DATA

rawFileName='F150818-0001E006.mat';
rawChannel='CRAW_006';
spkChannel='CSPK_006';
txtFileSortedLinear='F150818-0001E006.matCRAW_006FiltFilt.txt';

load (rawFileName);

rawSig=eval(rawChannel);
rawSig=double((rawSig*38.147)/200);
spkSig=eval(spkChannel);
spkSig=double((spkSig*38.147)/200);


HP=300;
LP=6000;
Fs=44000;
[b, a] = butter(4, [HP LP]/(Fs/2), 'bandpass');
nonlinearFilteredSig=filter(b, a, rawSig);

HPspk=300;
LPspk=6000;

[bAO, aAO] = butter(4, [HPspk LPspk]/(Fs/2), 'bandpass');

correctedAOspk=fliplr(filter(bAO,aAO,fliplr(spkSig)));


subplot(2,3,[4:5]);hold on;
plot ([0:length(spkSig)-1]/Fs,spkSig,'b');
plot ([0:length(rawSig)-1]/Fs,rawSig,'k');
plot ([0:length(correctedAOspk)-1]/Fs,correctedAOspk,'g');

axis ([1.55 1.75 -400 400]);
xlabel('Time (s)');
ylabel ('V (µV)');
title ('Electrophysiological data');

%WAVEFORM
x=csvread(txtFileSortedLinear,1);
linearWFtimeStamps=x(find(x(:,2)==1),3);
linearWFtimeStamps=round((linearWFtimeStamps)*44000);
timeWinBefore=1*44;
timeWinAfter=1.5*44;
rawWF=zeros(length(linearWFtimeStamps),timeWinBefore+timeWinAfter+1);
nonlinearWF=zeros(length(linearWFtimeStamps),timeWinBefore+timeWinAfter+1);
correctedWF=zeros(length(linearWFtimeStamps),timeWinBefore+timeWinAfter+1);

for i=1:length (linearWFtimeStamps)
    rawWF(i,:)=rawSig((linearWFtimeStamps(i)-timeWinBefore):(linearWFtimeStamps(i)+timeWinAfter));
    correctedWF(i,:)=correctedAOspk((linearWFtimeStamps(i)-timeWinBefore):(linearWFtimeStamps(i)+timeWinAfter));
    nonlinearWF(i,:)=nonlinearFilteredSig((linearWFtimeStamps(i)-timeWinBefore):(linearWFtimeStamps(i)+timeWinAfter));
end

subplot (2,3,6); hold on;
meanRawWF=(mean(rawWF))-10;
meanNonLinearWF=mean(nonlinearWF);
meancorrectedWF=mean(correctedWF)-20;
plot([0:size(meanRawWF,2)-1]/44,meanRawWF,'k');
plot([0:size(meanNonLinearWF,2)-1]/44,(meanNonLinearWF),'b');
plot([0:size(meancorrectedWF,2)-1]/44,(meancorrectedWF),'g');

title ('Waveform');
xlabel('Time (ms)');
ylabel ('V (µV)');






