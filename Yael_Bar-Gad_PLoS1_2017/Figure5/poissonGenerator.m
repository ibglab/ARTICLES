function [probVec]= poissonSpikeGenerator (T, rate, deltaT, refPer)
%T, deltaT, refPer [msec]
%rate [spikes/sec]

probVec=rand(T/deltaT,1);

i=1;
while i<=length(probVec)
    if probVec(i)<(rate)*(deltaT)*(1/1000)
        probVec(i)=1;
        if i+refPer<=length(probVec)
            probVec(i+1:i+refPer)=0;
            i=i+refPer;
        else
            probVec(i:end)=0;
            i=length(probVec);
        end
    else
        probVec(i)=0;
    end
    i=i+1;
end

