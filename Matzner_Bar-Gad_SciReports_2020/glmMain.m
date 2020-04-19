function [stimFilter, postSpikeFilter, rateBias, psthCorr] = glmMain(spikeTrain, stim, trainPercent, postSpikeFilterLength, stimFilterLength, baseBinWidth)
% Inputs:
% spikeTrain - n Trials x k bins
% stimlus - length k
% The fs of both inputs is 1000

% Defaults
if nargin<3    
    trainPercent = 0.8;
end
if nargin<4
    postSpikeFilterLength = 100;     % in ms
end
if nargin<5    
    stimFilterLength =60;       % in ms
end
if nargin<6    
    baseBinWidth = 4;   %boxcars bases
end
filtersSmoothWin = 9;
validateInitTime = 200;

T = length(stim);
trainTotalTime = floor(T*trainPercent);
trainStim = stim(1:trainTotalTime);
trainTrials = spikeTrain(:, 1:trainTotalTime);
testTrials = spikeTrain(:, trainTotalTime+1:end);

% Fit GLM
[rateBias, stimFilter, postSpikeFilter, trainPredictors, fval, weights, aic] = runGlm(trainStim, trainTrials, stimFilterLength, postSpikeFilterLength, filtersSmoothWin, baseBinWidth);
% Plot model parameters
plotGlmParams(stimFilter, postSpikeFilter, rateBias);
% Generate simulated spike trains and compare to real ones
[psthCorr, testPsth, modelPsth, h] = validateGlmModel(stimFilter, postSpikeFilter, rateBias, stim, testTrials, 'Test', validateInitTime);



