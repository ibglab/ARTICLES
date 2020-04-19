function [predictors, spikeTrain] = generatePredictors(stim, trials, stimBases, postSpikeBases)
% This functions gets the stimFilterLength postSpikeFilterLength and
% creates a matrix with the relevant predictors, including a column of
% [ones] for the bias term, and the appropriate spike train, using bases 
% stim = the stimulus 
% trials - matrix of N spike trains in the form [nTrials, length(stim)]
% stimFilterLength - in ms
% postSpikeFilterLength - in ms
% Output:
% Predictors - Matrix with stimFilterLength+postSpikeFilterLength+1 regressors X N
% spikeTrain = trials shaped to one vector of length N (nTrials*stimLength -minus
%                               beginning and end of trials, which are NaNs in the predictors matrix)

[nTrials, stimLength] = size(trials);
spikeTrain = reshape(trials', [nTrials*stimLength, 1]);

[stimFilterLength, nStimBases] = size(stimBases);
stimBinWidth = stimFilterLength/nStimBases;
[postSpikeFilterLength, nPostSpikeBases] = size(postSpikeBases);
postSpikeBinWidth = postSpikeFilterLength/nPostSpikeBases;
predictors = zeros(stimLength, nStimBases+nPostSpikeBases);
for i = 1:nStimBases
    s=conv(stim, stimBases(:, nStimBases-i+1));
    s(1:stimBinWidth*(nStimBases-i)) = NaN;
    s = [NaN, s(1:stimLength-1)];
    predictors(:, i) =s;
end
predictors = repmat(predictors, nTrials, 1);

for i = 1:nTrials
    trialRows = stimLength*(i-1);
    for j = 1:nPostSpikeBases        
        s=conv(trials(i,:), postSpikeBases(:, nPostSpikeBases - j + 1));
        s(1:postSpikeBinWidth*(nPostSpikeBases-j)) = NaN;
        s = [NaN, s(1:stimLength-1)];
        predictors([1: stimLength] + trialRows, j+nStimBases) =s;
    end
end

predictors = [ones(size(predictors, 1), 1), predictors];

% remove NaN values
noNanIdx = ~any(isnan(predictors), 2);
predictors = predictors(noNanIdx, :);
spikeTrain = spikeTrain(noNanIdx);

end

