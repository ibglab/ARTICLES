function []= figure3AB ()

Fs=44000;
T=120; %s

%waveform
load('STR_AP_3.mat');
samplingRateIC=20000;
ICspk=c011_AxoPatch(3650:3800);
derivativeIC=-diff(ICspk);
spkWF=interp (derivativeIC, round(Fs*10/samplingRateIC));
spkWF=resample(spkWF,1,10);
spkWF=spkWF*7;

%resample to lower width
spkWF=resample(spkWF,4,10);


[minWF,minWFindex]=min(spkWF);
timeWinBefore=(0.5*Fs/1000);
timeWinAfter=(1.5*Fs/1000);

spkWF=spkWF((minWFindex-timeWinBefore):(minWFindex+timeWinAfter)-1);
figure;
plot(spkWF);

[minWF,minWFindex]=min(spkWF);


% width
spkMat=zeros (5, length(spkWF));
counter=0;
for i=1:0.3:4.5
    counter=counter+1;
    spkWF2=interp (spkWF, i*10);
    spkWF2=resample(spkWF2,1,10);
    [minWF,minWFindex]=min(spkWF2);
    spkWF2=spkWF2(minWFindex-timeWinBefore:minWFindex+timeWinAfter-1);
    spkMat(counter,:)=spkWF2;  
    dftWF(counter,:)=fft(spkWF2);
end

%plot WFs
figure; hold on;
c = colormap((jet(size(spkMat,1)+10)));
c([1:5],:)=[];
c=flipud(c);
c([1:5],:)=[];
for ii = 1:(size(spkMat,1))
    hold on
    plot([0:length(spkWF)-1]/(Fs/1000),spkMat(ii,:),'Color',c(ii,:)');
end
hold off;
colorbar;



%plot dft
figure; hold on;
c = colormap((jet(size(dftWF,1)+10)));
c([1:5],:)=[];
c=flipud(c);
c([1:5],:)=[];
for ii = 1:(size(dftWF,1))
    hold on
    plot(abs(dftWF(ii,:)),'Color',c(ii,:)');
end
hold off;
colorbar;
v=axis;
axis ([2 12 v(3) v(4)]);

%spkTime
spkTrain=poissonGenerator (T*1000, 15, 1, 1);
spkTime=find(spkTrain==1)./1000; %s

% noise
pinkNoise=pinknoise (T*Fs);
pinkNoise=(pinkNoise./10)./2;

% LINEAR FILTER USING FILTFILT
%band-pass
HP =300;
LP =  6000;
[bbp, abp] = butter(4, [HP LP]/(Fs/2), 'bandpass');

% NON-LINEAR
[b, a] = butter(4, [HP LP]/(Fs/2), 'bandpass');



for j=1:size(spkMat,1)
    
    timeWinAfter=(1*Fs/1000);
    
    %spk signal
    spkSig=zeros(1,T*Fs);
    x=round(spkTime*Fs);
    for i=1:length (x)
        spkSig ((x(i)-timeWinBefore):(x(i)-timeWinBefore+length(spkWF)-1))=spkSig ((x(i)-timeWinBefore):(x(i)-timeWinBefore+length(spkWF)-1))+spkMat(j,:);
    end
    
    %spk+noise
    signal=pinkNoise+spkSig;
    
    %linear filter
    linearFilteredSig=filtfilt(bbp, abp, signal);
    
    %non linear filter
    nonlinearFilteredSig=filter(b, a, signal);
    
    
    %meanWF
    rawWF=[];
    linearWF=[];
    nonlinearWF=[];
    for i=1:length (spkTime)
        rawWF(i,:)=signal((x(i)-timeWinBefore):(x(i)+timeWinAfter-1));
        linearWF(i,:)=linearFilteredSig((x(i)-timeWinBefore):(x(i)+timeWinAfter-1));
        nonlinearWF(i,:)=nonlinearFilteredSig((x(i)-timeWinBefore):(x(i)+timeWinAfter-1));
    end
    meanRawWF=(mean(rawWF));
    meanRawWF=meanRawWF-meanRawWF(1);
    meanNonLinearWF=mean(nonlinearWF);
    meanLinearWF=(mean(linearWF));
    meanLinearWF=meanLinearWF-meanLinearWF(1);
    
    %recalc aligned WF
    
    [vRaw,iRaw]=min(meanRawWF);
    [vNLP,iNLP]=min(meanNonLinearWF);
    [vLP,iLP]=min(meanLinearWF);
    
    diffNLPWF=iNLP-iRaw;
    %diffLPWF=iLP-iRaw;
    diffLPWF=0;
    
    linearWF=[];
    nonlinearWF=[];
    for jj=1:length (spkTime)
        linearWF(jj,:)=linearFilteredSig((x(jj)+diffLPWF-timeWinBefore):(x(jj)+diffLPWF+timeWinAfter-1));
        nonlinearWF(jj,:)=nonlinearFilteredSig((x(jj)+diffNLPWF-timeWinBefore):(x(jj)+diffNLPWF+timeWinAfter-1));
    end
    
    
    meanNonLinearWF=mean(nonlinearWF);
    meanLinearWF=(mean(linearWF));
    meanLinearWF=meanLinearWF-meanLinearWF(1);
    
    
    figure; hold on;
    plot([0:size(meanRawWF,2)-1]/(Fs/1000),meanRawWF,'k');
    plot([0:size(meanLinearWF,2)-1]/(Fs/1000),(meanLinearWF),'r');
    plot([0:size(meanNonLinearWF,2)-1]/(Fs/1000),(meanNonLinearWF),'b');
    DNonLinear(j)=sqrt(sum ((meanNonLinearWF-meanRawWF).^2))/sqrt( sum ((meanRawWF).^2));
    DLinear(j)=sqrt(sum ((meanLinearWF-meanRawWF).^2))/ sqrt(sum ((meanRawWF).^2));

end

figure; hold on;
plot (1:j, DNonLinear,'b*-');
plot (1:j, DLinear,'r*-');
axis ([1 j 0.1 1.3]);

