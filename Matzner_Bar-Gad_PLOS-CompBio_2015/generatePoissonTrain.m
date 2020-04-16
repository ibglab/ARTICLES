function [spikeTrain] = generatePoissonTrain(totalTime, rateFunction, refractoryPeriod)
%Generate poisson spike train using a rate function

if (refractoryPeriod == 0)
    spikeTrain = rand(1,totalTime) < rateFunction(1:totalTime);
else
     spikeTrain=zeros(1,totalTime);
    for i=[1:totalTime]
        if (rand<rateFunction(i))
        spikeTrain(i)=1;
            for j=[1:refractoryPeriod]       % set refactory period rate to 0
              rateFunction(i+j)= 0;
            end       
        end
    end
end
end