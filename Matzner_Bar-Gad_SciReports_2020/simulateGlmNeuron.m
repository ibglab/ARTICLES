function [spikeTrain, Vstim, VpostSpike, Vmem] = simulateGlmNeuron(stim, stimFilter, postSpikeFilter, rateBias)
% Simulates a GLM neuron, where the post spike filter is a modulation of
% the base firing rate
% stim  - the stimulus
%   stimFilter - the stimulus filter
%   postSpikeFilter = the post spike filter
% rateBias - the bias term (in the format of log(spikes/s / 1000)
% Output:
% spikeTrain
% Vstim - voltage resulted from stimulus
% VpostSpike - voltage resulted from spiking activity
% Vmem - the total voltage

startTime = 1;      % start simulation from the beginning of stimulus  (the beginning of the spike train will not contain the full model stim filters)

if (size(stim, 1) ==1)
    stim = stim';
end
postSpikeFilterLength = length(postSpikeFilter);
totalTime = length(stim);

spikeTrain = zeros(totalTime, 1);
rate = zeros(totalTime, 1);
VpostSpike = zeros(totalTime,1);
Vstim = conv(stim ,fliplr(stimFilter));
Vstim = [0;Vstim(1:totalTime-1)] + rateBias;

for t =startTime:totalTime
    rate(t) = exp(Vstim(t)+ VpostSpike(t));
    if rand < rate(t)
        spikeTrain(t) = 1;
        postSpikeBins = min(postSpikeFilterLength, totalTime-t);
        %              VpostSpike(t+(1  :postSpikeBins)) =  postSpikeFilter(1:postSpikeBins)';              %Reset the neuron
        VpostSpike(t+(1:postSpikeBins)) = VpostSpike(t+(1:postSpikeBins)) + postSpikeFilter(1:postSpikeBins)';        %Adds the filter to the current state
    end
end

Vstim = Vstim - rateBias;
Vmem = log(rate);