function [psthCorr, testPsth, modelPsth, h] = validateGlmModel(stimFilter, postSpikeFilter, rateBias, stim, testTrials, cellName, initTime)
%   Gets GLM parameters, and simulates new spike trains. Then compares the
%   simulated PSTH to the test data PSTH using correlation coefficient.
%   initStim - time to initialize the simulation before strating the recording
%  stim - full stim (train + test)

if (nargin < 7)
    initTime = 0;
end
[nTrials, testLength] = size(testTrials);
stimLength = length(stim);
trainLength = stimLength - testLength;
testTimes = [trainLength+1 : stimLength];

h=figure; 
ax1= subplot(4, 1, 1); 
plot(testTimes, stim(testTimes));  set(gca, 'xtick', []); ylabel('Stim');

ax2 = subplot(4, 1, 2:3);
% Plot test rasters
for i = 1:nTrials
    spikeTimes = find(testTrials(i, :));
    plot(spikeTimes+trainLength, ((nTrials*2)-(i-1))*ones(1, length(spikeTimes)), 'k.'); hold on;
end

% Simulate and plot model rasters
stim = stim(end - testLength - initTime +1 : end);
modelTrials = zeros(nTrials, testLength);
for i = 1:nTrials
    [spikeTrain, Vstim, VpostSpike, Vmem] = simulateGlmNeuron(stim, stimFilter, postSpikeFilter, rateBias);
    spikeTrain = spikeTrain(initTime+1:end);
    modelTrials(i, : )= spikeTrain;
    spikeTimes = find(spikeTrain);
    plot(spikeTimes+trainLength, (nTrials-i+1)*ones(1, length(spikeTimes)), 'r.'); hold on;
end
ylim([0 nTrials*2]);
yLabels = get(gca, 'yticklabel');
set(gca, 'yticklabel', flipud(yLabels));  ylabel('Trials');
xlim([testTimes(1)-1 testTimes(end)]);

testPsth = smoothedPsth(testTrials, 11).*1000;
modelPsth = smoothedPsth(modelTrials, 11).*1000;
psthCorr = corr(testPsth', modelPsth');
if (isnan(psthCorr))        % Occurs mostly when vectors are both zeros (no spikes)
    psthCorr = 0;
end
psthMse = goodnessOfFit(testPsth, modelPsth, 'MSE');

ax3 = subplot(4, 1,4); 
plot(testTimes, testPsth, 'k'); hold on; plot(testTimes, modelPsth, 'r');
xlabel('Time (ms)'); ylabel('PSTH');
title(['Corr = ', num2str(psthCorr)]);
suptitle(cellName);

linkaxes([ax1,ax2,ax3], 'x');
