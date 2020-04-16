function [] = plotStimFilter(stimFilter)
% Plots stim filter properly - where the x-axis is Time from -length of
% filter until Time 0
% stimFilter- can be a single stim filter, or multiple filter in an array [lengthOfFilter Filters]
timeInterval = 5;       %Intervals of time in the x-axis
stimFilterLength = size(stimFilter, 2);
 plot([0:stimFilterLength-1], stimFilter);title('Stimulus Filter'); 
tickTimes = [stimFilterLength:-round(stimFilterLength/timeInterval):0]; 
set(gca, 'xtick', fliplr(tickTimes), 'xticklabel', fliplr(tickTimes - stimFilterLength));
xlabel('Time (ms)'); 
xlim([0 stimFilterLength]);

end

