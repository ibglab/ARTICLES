function [truePositive, falsePositive]= ROC(signal,spkTime, Fs)

threshold=[2:0.25:10];%stds
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




