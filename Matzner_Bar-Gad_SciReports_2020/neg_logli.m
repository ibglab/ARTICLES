function [neLogLi, grad, hess ] = neg_logli( predictors, weights, spikeTrain)
% copmutes negative log-likelihood for  of spike train with Poisson model
% and exponential non-linarity
% Inputs:
% predictos -  matrix of M regressors in the form [N x numPredictors]   (N is the total length of the stimulus)
% weights - regression weights [M x 1]
% spikeTrain - output spike train [N x 1]
E=predictors*weights; 
rate=exp(E);

zIdx = rate==0;
if (find(zIdx))
    warning('Zero values in rate function');
    rate(find(zIdx)) = exp(-745);       % the minimum before getting -Inf
end

neLogLi = -sum(E(find(spikeTrain))) +sum(rate);
grad = predictors' * (rate - spikeTrain);   % this grad is derived using the first option of  neLogLi
 hess = predictors'*bsxfun(@times,predictors,rate); % non-spiking term 

end

