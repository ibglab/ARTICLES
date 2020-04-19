function []= spikeSortingSimulatedData ()

threshold=2.5;%stds


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
pinkNoise=(pinkNoise./10)./2;

signal=pinkNoise+spkSig;

% LINEAR FILTER USING FILTFILT
%band-pass
HP =300;
LP =  6000;
[b, a] = butter(4, [HP LP]/(Fs/2), 'bandpass');
linearFilteredSig=filtfilt(b, a, signal);

% NON-LINEAR
nonlinearFilteredSig=filter(b, a, signal);

figure; hold on;
plot ([0:length(linearFilteredSig)-1]/Fs,linearFilteredSig);
plot ([[0,length(linearFilteredSig)-1]/Fs],[mean(linearFilteredSig)-threshold*std(linearFilteredSig),mean(linearFilteredSig)-threshold*std(linearFilteredSig)],'g');


% detection
[x,correctDetectionLinear]=detection (linearFilteredSig,threshold,Fs, spkTime);
[y,correctDetectionNonlinear]=detection (nonlinearFilteredSig,threshold,Fs,spkTime);

%meanWF linear
for i=1:length (x)-1
    linearWF(i,:)=linearFilteredSig((x(i)-30):(x(i)+80));
end

%meanWF nonlinear
for i=1:length (y)-1
    nonlinearWF(i,:)=nonlinearFilteredSig((y(i)-30):(y(i)+80));
end

pca (linearWF,correctDetectionLinear(1:end-1));
title ('LP filter');
pca (nonlinearWF,correctDetectionNonlinear(1:end-1));
title ('NLP filter');





function [x,correctDetection] = detection (signal, threshold,Fs, spkTime)

refPeriod=1*Fs/1000;%ms
WFwidth=0.001;%s

meanSignal=mean(signal);
stdSignal=std(signal);
thresholdSignal=meanSignal-threshold.*stdSignal;

y=[refPeriod:refPeriod:length(signal)-1];

for j=1:length(thresholdSignal)
    %detection
    x=find(signal<=thresholdSignal(j));
    detection=zeros(1,length(x));
    i=1;
    while i<length (x)
           detection(i)=1;
           i=find (x>(x(i)+refPeriod),1);      
    end
    x=x(detection==1);
    
    % find true positive (detection of real spikes)
    correctDetection=zeros(1,length(x));
    detectedSpkTimes=[];
    missDetectedSpk=zeros(1,length(x));
    for i=1:length(spkTime)
        a=find ((x>=(spkTime(i)-WFwidth)*Fs) & (x<(spkTime(i)+WFwidth)*Fs));
        if (length(a)>0)
            correctDetection(a)=1;
            detectedSpkTimes(i)=1;
        else
            missDetectedSpk(a)=1;
        end
    end
    detectedSpk=spkTime(detectedSpkTimes==1);
    missDetectedSpk=spkTime(find(missDetectedSpk==1));
    truePositive(j)= length(detectedSpk)/length(spkTime);
    
    
    % find false positive (detected non-spikes)
    falsePositive (j)=(length(x)-length(detectedSpk))/(length(y)-length(spkTime));
    
end



function []= pca (WF, truePositive)
[w,pc]=princomp(WF);
pc=pc';

X=([pc(1,:);pc(2,:)]');
[idx,C]=classification (X);

figure;
subplot(3,3,[4:9]);hold on;
plot (X(idx==1,1),X(idx==1,2),'r.');
plot (X(idx==2,1),X(idx==2,2),'b.');
plot (pc(1,find(truePositive==1)),pc(2,find(truePositive==1)),'ko');
plot (C(:,1),C(:,2),'kx');
xlabel('PC1');
ylabel('PC2');
legend ('Cluster1', 'Cluster2', 'Real spikes')

subplot(3,3,1:2); hold on;
plot (WF(idx==1,:)','r');
plot (WF(idx==2,:)','b');

plot (mean(WF(idx==1,:)),'k','Linewidth',2);
plot (mean(WF(idx==2,:)),'k','Linewidth',2);
plot (mean(WF(idx==1,:))+std(WF(idx==1,:)),'k:');
plot (mean(WF(idx==1,:))-std(WF(idx==1,:)),'k:');
plot (mean(WF(idx==2,:))+std(WF(idx==2,:)),'k:');
plot (mean(WF(idx==2,:))-std(WF(idx==2,:)),'k:');
axis tight;



function [idx,C]=classification (X)

[idx,C]=kmeans (X,2);



