function []= figure3CD ()
%cutoffs

Fs=44000;
T=120; %s

%waveform
load('STR_AP_3.mat');
samplingRateIC=20000;
ICspk=c011_AxoPatch(3650:3800);
derivativeIC=-diff(ICspk);
spkWF=interp (derivativeIC, round(Fs*10/samplingRateIC));
spkWF=resample(spkWF,1,10);
spkWF=spkWF*10;

[minWF,minWFindex]=min(spkWF);
timeWinBefore=(0.4*Fs/1000);
timeWinAfter=(0.5*Fs/1000);
spkWF=spkWF((minWFindex-timeWinBefore):(minWFindex+timeWinAfter)-1);
figure;
plot(spkWF);


%spkSig
spkTrain=poissonGenerator (T*1000, 15, 1, 1);
spkTime=find(spkTrain==1)./1000; %s
x=round(spkTime*Fs);
spkSig=zeros(1,T*Fs);
for i=1:length (x)
    spkSig ((x(i)-timeWinBefore):(x(i)-timeWinBefore+length(spkWF)-1))=spkWF;
end

% noise
pinkNoise=pinknoise (T*Fs);
pinkNoise=(pinkNoise./10)./2;

%spk+noise
signal=pinkNoise+spkSig;



 
% LINEAR FILTER USING FILTFILT
%band-pass
HP =[50:50:500];
LP =  6000;

for i=1:length(HP)
[b, a] = butter(4, [HP(i) LP]/(Fs/2), 'bandpass');

%LP filter
linearFilteredSig=filtfilt(b, a, signal);

%NLP filter
nonlinearFilteredSig=filter(b, a, signal);


    
%meanWF
rawWF=zeros(length(spkTime),length(spkWF));
linearWF=zeros(length(spkTime),length(spkWF));
nonlinearWF=zeros(length(spkTime),length(spkWF));
for jj=1:length (spkTime)
    rawWF(jj,:)=signal((x(jj)-timeWinBefore):(x(jj)+timeWinAfter-1));
    linearWF(jj,:)=linearFilteredSig((x(jj)-timeWinBefore):(x(jj)+timeWinAfter-1));
    nonlinearWF(jj,:)=nonlinearFilteredSig((x(jj)-timeWinBefore):(x(jj)+timeWinAfter-1));
end
meanRawWF=(mean(rawWF));
meanNLPWF=mean(nonlinearWF);
meanLPWF=(mean(linearWF));
meanLPWF=meanLPWF-meanLPWF(5);

%align WF
[vRaw,iRaw]=min(meanRawWF);
[vNLP,iNLP]=min(meanNLPWF);
[vLP,iLP]=min(meanLPWF);

diffNLPWF=iNLP-iRaw;
diffLPWF=0;

linearWF=zeros(length(spkTime),length(spkWF));
nonlinearWF=zeros(length(spkTime),length(spkWF));
for jj=1:length (spkTime)
    linearWF(jj,:)=linearFilteredSig((x(jj)+diffLPWF-timeWinBefore):(x(jj)+diffLPWF+timeWinAfter-1));
    nonlinearWF(jj,:)=nonlinearFilteredSig((x(jj)+diffNLPWF-timeWinBefore):(x(jj)+diffNLPWF+timeWinAfter-1));
end

meanRawWF(i,:)=(mean(rawWF));
meanNonLinearWF(i,:)=mean(nonlinearWF);
meanLinearWF(i,:)=(mean(linearWF));
meanLinearWF(i,:)=meanLinearWF(i,:)-meanLinearWF(i,5);

figure;hold on;
plot (mean(rawWF),'k');
plot (mean(nonlinearWF),'b');
plot (meanLinearWF(i,:)-meanLinearWF(i,1),'r');
title (['cutoff=' num2str(HP(i)) 'Hz'])


DNonLinear(i)=sqrt(sum ((meanNonLinearWF(i,:)-meanRawWF(i,:)).^2))/sqrt( sum ((meanRawWF(i,:)).^2));
DLinear(i)=sqrt(sum ((meanLinearWF(i,:)-meanRawWF(i,:)).^2))/ sqrt(sum ((meanRawWF(i,:)).^2));
end


figure; hold on;
plot([0:size(meanRawWF,2)-1]/(Fs/1000),meanRawWF,'k');


c = colormap(flipud(autumn(size(meanLinearWF,1))));
for ii = 1:(size(meanLinearWF,1))
    plot([0:size(meanLinearWF,2)-1]/(Fs/1000),(meanLinearWF(ii,:)),'Color',c(ii,:)');
end

c = colormap(flipud(winter(size(meanNonLinearWF,1)+5)));
c([1:5],:)=[];
for ii = 1:(size(meanNonLinearWF,1))
    plot([0:size(meanNonLinearWF,2)-1]/(Fs/1000),(meanNonLinearWF(ii,:)),'Color',c(ii,:)');
end
hold off;

figure;
cLP = colormap(flipud(autumn(size(meanLinearWF,1))));
colorbar;

figure;
cNLP = colormap(flipud(winter(size(meanNonLinearWF,1)+5)));
cNLP([1:5],:)=[];
colorbar;

figure; hold on;
plot (1:i, DNonLinear,'b*-');
plot (1:i, DLinear,'r*-');


figure; hold on;
colormapline (1:i, DNonLinear,[],cNLP);
colormapline (1:i, DLinear,[], cLP);
