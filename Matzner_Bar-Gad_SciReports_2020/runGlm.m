function [ rateBias, stimFilter, postSpikeFilter, trainPredictors, fval, b, aic ] = runGlm(stim, trials, stimFilterLength, postSpikeFilterLength, smoothWin, baseBinWidth)
% stim = the stimulus 
% trials - matrix of N spike trains in the form [nTrials, length(stim)]
% stimFilterLength - in ms
% postSpikeFilterLength - in ms
% smoothWin - win size for smoothing filters (must be odd!). For no smoothing insert 1.
% baseBinWidth - 1 indictes no base
 
if (nargin < 6)
    baseBinWidth = 1;
end

if (baseBinWidth > 1)
    [stimBases, nStimBases] = generateEqualBoxcarBases(stimFilterLength, baseBinWidth);
    [postSpikeBases, nPostSpikeBases]  = generateEqualBoxcarBases(postSpikeFilterLength, baseBinWidth);
    [trainPredictors, trainSpikeTrain] = generatePredictors(stim, trials, stimBases, postSpikeBases);
else 
    [trainPredictors, trainSpikeTrain] = generatePredictorsWithoutBases(stim, trials, stimFilterLength, postSpikeFilterLength);
end

fun = @(w)(neg_logli(trainPredictors, w, trainSpikeTrain));
w0 = trainPredictors\trainSpikeTrain; 
opts = optimset( 'GradObj', 'on', 'Hessian', 'on', 'display', 'iter');
tic
[b,fval,exitFlag, output, grad, hess]= fminunc(fun, w0, opts);
toc

rateBias = b(1);
if (baseBinWidth > 1)
    stimFilter = stimBases * b(2:nStimBases+1);
    postSpikeFilter = flipud(postSpikeBases * b(nStimBases+2:end));
else
    stimFilter = b(2:stimFilterLength+1)';
    postSpikeFilter = fliplr(b(stimFilterLength+2:end)');
end

% Smooth filters
stimFilter = yySmoothData(stimFilter, 'gauss', smoothWin);
postSpikeFilter= yySmoothData(postSpikeFilter, 'gauss', smoothWin);

% Akai information criterion
aic = 2*fval + 2*length(b);

end

