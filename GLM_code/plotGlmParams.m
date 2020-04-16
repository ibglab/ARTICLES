function h = plotGlmParams( stimFilter, postSpikeFilter, rateBias, cellName)
% Plots stim filter, post spike filter(exponential) and rate bias of a single cell
postSpikeFilterLength = length(postSpikeFilter);
h = figure; 
subplot(1, 3, 1);
plotStimFilter(stimFilter);
subplot(1, 3, 2); plot(exp(postSpikeFilter)); title('Post Spike Filter'); 
tickTimes = [0:round(postSpikeFilterLength/5):postSpikeFilterLength]; 
set(gca, 'xtick', tickTimes); xlabel('Time (ms)'); xlim([0 postSpikeFilterLength]);
bias = exp(rateBias)*1000;
subplot(1, 3, 3); plot(1, bias, '.'); title('Bias'); ylim([0 bias+10]);
set(gca, 'xticklabel', []);
if (nargin > 3)
    suptitle(cellName);
end

