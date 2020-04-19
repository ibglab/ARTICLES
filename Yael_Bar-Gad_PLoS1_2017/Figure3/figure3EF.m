function []= figure3EF ()

Fs=44000;
T=120; %s

%waveform
load('D:\Dorin\filters paper\IC data\STR_AP_3.mat');
samplingRateIC=20000;
ICspk=c011_AxoPatch(3650:3800);
derivativeIC=-diff(ICspk);
spkWF=interp (derivativeIC, round(Fs*10/samplingRateIC));
spkWF=resample(spkWF,1,10);
spkWF=spkWF*7;

%spkTime
spkTrain=poissonGenerator (T*1000, 15, 1, 1);
spkTime=find(spkTrain==1)./1000; %s

% noise
pinkNoise=pinknoise (T*Fs);
pinkNoise=(pinkNoise./10)./2;
stdPinkNoise=std(pinkNoise);


sigAmp=[0.2:0.2:2];
spkTemp=spkWF;
counter=0;
for kk=sigAmp
    counter=counter+1;
    spkWF=spkTemp*kk;
    snr(counter)= max(abs(spkWF))./stdPinkNoise;
    
    %spk signal
    spkSig=zeros(1,T*Fs);
    x=round(spkTime*Fs);
    for i=1:length (x)
        spkSig ((x(i)-80):(x(i)-80+length(spkWF)-1))=spkSig ((x(i)-80):(x(i)-80+length(spkWF)-1))+spkWF;
    end
    
    signal=pinkNoise+spkSig;
       
    % LP FILTER
    %band-pass
    HP =300;
    LP =  6000;
    [bbp, abp] = butter(4, [HP LP]/(Fs/2), 'bandpass');
    linearFilteredSig=filtfilt(bbp, abp, signal);
    meanLinear=mean(linearFilteredSig);
    stdLinear=std(linearFilteredSig);
    
    % NLP FILTER
    [b, a] = butter(4, [HP LP]/(Fs/2), 'bandpass');
    nonlinearFilteredSig=filter(b, a, signal);
    
    %ROC
    [truePositiveNon, falsePositiveNon]= ROC(nonlinearFilteredSig,spkTime, Fs);
    [truePositive, falsePositive]= ROC(linearFilteredSig,spkTime, Fs);
    
    truePositive=[fliplr(truePositive), 1];
    falsePositive=[fliplr(falsePositive),1];
    truePositiveNon=[fliplr(truePositiveNon), 1];
    falsePositiveNon=[fliplr(falsePositiveNon),1];
    
    %plot ROC 
    figure; hold on;
    plot (falsePositive,truePositive,'r');
    plot (falsePositiveNon,truePositiveNon);
    plot ([0,1],[0,1],'k');
    ylabel ('P(1-ß)');
    xlabel ('P(a)');
    title (['SNR=' num2str(snr(counter))]);
    
    %calc area under the ROC curve
    Q(counter)=trapz((falsePositive),(truePositive));
    Qnon(counter)=trapz((falsePositiveNon),(truePositiveNon));
      
end

figure;hold on;
plot (snr,Q,'r');
plot (snr,Qnon);
xlabel('Signal to Noise Ratio');
ylabel ('? (1-ß)da');
xlim ([0.5 5.5]);
