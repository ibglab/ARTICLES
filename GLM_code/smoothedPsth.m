function [ psth ] = smoothedPsth( trials , winSize)
% Computes smoothed rate of spike trains during stimulus over N trials.
% trials - matrix of N spike trains in the form [nTrials, length(stim)]
spMatAll = sum(trials,1);
psth = yySmoothData(spMatAll,'gauss',winSize)/size(trials,1);

end

