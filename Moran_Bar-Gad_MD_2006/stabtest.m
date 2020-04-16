function res = stabtest(sig, len, win, alpha)
% res = stabtest(sig, len, win, alpha)
%
% Checks the stability of an analog signal using ANOVA test on RMS values 
% of windows on start and end of signal. 
% sig - the analog signal, in microvolts.
% len - the length of samples to chaeck on  both sides of the signal
% win - the window size in samples. len should be multiples of win
% alpha - the rejection limit

% in the article the followitn parameters were used:
% sig - 5 seconds analog signal samples at 24KHz
% len - 2 soconds * 24KHz samples
% win - 20 miliseconds * 24KHz / 1000
% alpha - 0.005


    % there should be at least 2 distinct samples to evaluate
    if (length(sig) < 2*len)
        res = -1;
        return
    else
        % the first part windowed
        fst = reshape(sig(1:len), win, len/win);
        % the last part windowed
        lst = reshape(sig(end-len+1:end),win, len/win);
        
        % calculate the RMS of the windows in the two part
        fstrms = sqrt(sum(fst.^2)/win); 
        lstrms = sqrt(sum(lst.^2)/win);
        
        % evaluate the difference between the population means 
        P = anova1([fstrms' lstrms'], [], 'off');
        if P < alpha
            res = 1;
        else
            res = 0;
        end
    end