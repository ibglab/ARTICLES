function []= figure5D ()

Fs=44000;
T=120; %s

[spkTime, pinkNoise, spkWF]= simulation(Fs, T);


% calculate spike WFs (spkMat)
width=[1.6:0.4:4.5];

[minWF,minWFindex]=min(spkWF);
timeWinBefore=(0.8*Fs/1000);
timeWinAfter=(1*Fs/1000);

spkWF=spkWF((minWFindex-timeWinBefore):(minWFindex+timeWinAfter-1));
spkMat=zeros (length(width), length(spkWF));
counter=0;
for i=width
    counter=counter+1;
    spkWF2=interp (spkWF, round(i*10));
    spkWF2=resample(spkWF2,1,10);
    [minWF,minWFindex]=min(spkWF2);
    spkWF2=spkWF2(minWFindex-timeWinBefore:minWFindex+timeWinAfter-1);
    spkMat(counter,:)=spkWF2;
end

for i=1:size(spkMat,1) %runs over spike widths
    %calc signals
    spkWF=spkMat(i,:);
    [rawSig, NLPsig]= getSig(spkTime, pinkNoise, spkWF, Fs, T,timeWinBefore);
    
    %segment signal (segment size- maxSegment)
    maxSegment=round(2*Fs/1000);%ms
    minSegment=round(1*Fs/1000);%ms
    [rawWF, NLPWF] =segmentSig (rawSig, NLPsig, spkTime, maxSegment,minSegment, Fs);
    
    counter=0;
    bin=round(0.1*Fs/1000);
    while (size(NLPWF,2)>minSegment) %runs over segment width (shortens it every enterance)
        l=size(NLPWF,2);
        rawWF(:,l-bin:l)=[];
        NLPWF(:,l-bin:l)=[];
        counter=counter+1;
        
        %correct WF
        [correctedWF]= correctSegment(NLPWF,Fs,length(spkTime));
        
        %calc mean WFs
        meanRawWF=mean(rawWF);
        meanNLPWF=mean(NLPWF);
        meanCorrectedWF=mean(correctedWF);
        meanCorrectedWF=meanCorrectedWF-meanCorrectedWF(1)+meanRawWF(1);
        %calc D
        D(i,counter)=sqrt(sum ((meanCorrectedWF(1:minSegment-3)-meanRawWF(1:minSegment-3)).^2))/sqrt( sum ((meanRawWF(1:minSegment-3)).^2));
        
        %plot WFs
        figure; hold on;
        plot([0:size(meanRawWF,2)-1]/(Fs/1000),meanRawWF,'k');
        plot([0:size(meanNLPWF,2)-1]/(Fs/1000),(meanNLPWF),'b');
        plot([0:size(meanCorrectedWF,2)-1]/(Fs/1000),(meanCorrectedWF),'g');
        title (['WF width=' num2str(width(i)) ', seg width=' num2str(size(meanCorrectedWF,2)) ', D=' num2str(D(i,counter))])
        v=axis;
        axis ([0 2.5 v(3) v(4)]);
    end
end

D=fliplr(D);
figure;
plot(D);
ylabel('D');
xlabel('WF width');
ylim([0 0.25]);


%figure;hold on;
c = colormap((jet(size(D,1)+11)));
c([1:5],:)=[];
c=flipud(c);
c([1:5],:)=[];
for ii = 1:(size(D,1))
    hold on
    plot((D(:,ii)),'Color',c(ii,:)');
end
colorbar;
ylabel('D');
xlabel('WF width');
ylim([0 0.25]);
a=1;


end


function [spkTime, pinkNoise, spkWF]= simulation (Fs, T)
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

%spkTime
spkTrain=poissonGenerator (T*1000, 15, 1, 1);
spkTime=find(spkTrain==1)./1000; %s

% noise
pinkNoise=pinknoise (T*Fs);
pinkNoise=(pinkNoise./10)./2;
end

function [rawWF, NLPWF] =segmentSig (rawSig, NLPsig, spkTime, maxSegment,minSegment, Fs)


    timeWinBefore=0.4*Fs/1000;%ms
    timeWinAfter=maxSegment-timeWinBefore+1;
    
    %meanWF
    rawWF=[];
    NLPWF=[];
    x=round(spkTime*Fs);
    for i=2:length (spkTime)-1
        rawWF(i,:)=rawSig((x(i)-timeWinBefore):(x(i)+timeWinAfter-1));
        NLPWF(i,:)=NLPsig((x(i)-timeWinBefore):(x(i)+timeWinAfter-1));
    end

end



function [rawSig, NLPsig]= getSig(spkTime, noise, spkWF, Fs, T,timeWinBefore)

%spk signal
spkSig=zeros(1,T*Fs);
x=round(spkTime*Fs);
for i=1:length (x)
    spkSig ((x(i)-timeWinBefore):(x(i)-timeWinBefore+length(spkWF)-1))=spkSig ((x(i)-timeWinBefore):(x(i)-timeWinBefore+length(spkWF)-1))+spkWF;
end

%spk+noise
rawSig=noise+spkSig;

% NLP filter
HP =300;
LP =  6000;
[b, a] = butter(4, [HP LP]/(Fs/2), 'bandpass');
NLPsig=filter(b, a, rawSig);
end


function [correctedWF]= correctSegment(NLPWF,Fs, spkNum)

opt = 1;
l=size(NLPWF,2);
for k=1:size(NLPWF,1)
    switch (opt)
        case 1 %linear decay
            sl = (NLPWF(k,l)-NLPWF(k,l-2))/2;
            NLPWF(k,l+[1:10]) = NLPWF(k,l)+sl*[1:10];
        case 2 %constant value
            NLPWF(k,l+[1:10]) = NLPWF(k,l);
        case 3 %zero padding
            NLPWF(k,l+[1:10]) = 0;
        case 4%spline fit
            NLPWF(k,l+[1:10])  = interp1([1:l l+11:l+20],[NLPWF(k,1:l) zeros(1,10)],[l+1:l+10],'spline');
    end
end

%correct concatenated segments
concatSeg=reshape (NLPWF', ([1, size(NLPWF,1)*size(NLPWF,2)]));
correctedSig= correction (concatSeg, Fs);
[correctedWF] = (reshape(correctedSig,[size(NLPWF,2),size(NLPWF,1)]))';

%elongation (padding) removal
correctedWF(:,(l+1:end))=[];

end

function [correctedSig]= correction (sig, Fs)

HP=300;
LP=6000;
[b, a] = butter(4, [HP LP]/(Fs/2), 'bandpass');

correctedSig=filter(b,a, fliplr(sig));
correctedSig= fliplr(correctedSig);
end







