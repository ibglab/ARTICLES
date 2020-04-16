function [stnin stnout] = estimLocRMS(EDT, RMS, posterior, printit)
% [stnin stnout] = estimLocRMS(EDT, RMS, posterior, printit)
% This function estimate the entry and exit points of the STN based on
% liklihood probabilities calculated from previous surgeries

% Input parameters:
% EDT: array of estimated distance to target values (depth)
% RMS: array of normalized RMS values
% posterior: a structure containing the 3 2D posterior probabilities for each of the 
% locations - pre, within and post STN
% printit: a flag indicating wheather to print the graph of probabilities.
% 1 to print, 0 not to print

    % The 2D posterior probabilities for each of the locations
    p_stn_given_edt_rms = posterior.p_stn_given_dep_rms ;
    p_pre_given_edt_rms = posterior.p_pre_given_dep_rms ;
    p_post_given_edt_rms = posterior.p_post_given_dep_rms ;
    
    % the bins of the 2D plane of EDT and RMS
    rmsbins = posterior.rmsbins;
    depthbins = posterior.depthbins ;
    
    
    pstn = [];
    ppre = [] ;
    ppost = [] ;
    rej = [];
    count = 0 ;
    
    % loop through all EDT-RMS pairs and find their corresponding posterior
    % probabilities in each of the locations
    for i = 1:length(EDT)
        rmsbin = find(histc(RMS(i),rmsbins));
        depbin = find(histc(EDT(i), depthbins)) ;
        % reject sessions which are out of the valid paremeters range 
        if isempty(rmsbin) | isempty(depbin)
            rej = [rej i];
            continue
        end
        count = count + 1 ;
        % find the posterior probability being in each of the location
        % given the RMS and EDT
        pstn(count) = p_stn_given_edt_rms(rmsbin, depbin);
        ppre(count) = p_pre_given_edt_rms(rmsbin, depbin);
        ppost(count) = p_post_given_edt_rms(rmsbin, depbin);    
    end
    % remove the rejected sessions
    RMS(rej) = [] ;
    EDT(rej) = [] ;
    
    % sort the session in descending order
    [depsort,depind] = sort(EDT, 'descend');
    
    % find the places where the probability to be within the STN is higher
    % than being dorsal to it (preior)
    optin = find(ppre(depind) < pstn(depind));
    % find the first places where there are 2 succesive decisions of being 
    % within the STN
    ent = find(diff(optin)==1);
    stnin = EDT(optin(ent(1))); 
    
    % take the last place where the probability to be within the STN is 
    % still higher than being after (ventral) to the STN
    optout = find(pstn(depind) > ppost(depind));
    if (optout(end) == length(EDT))
        ind = optout(end);
    else
        ind = optout(end)+1;
    end
    stnout = EDT(ind);

    % print the RMS and the system's borders decisions
    if printit
        figure
        ax1 = gca ;
        h1 = bar(ax1, depsort,RMS(depind));
        set(h1, 'FaceColor', 'w')
        set(ax1, 'XDir', 'reverse'); %, 'YLim', [0.4e-5 2.8e-5]
        set(ax1,'XAxisLocation','bottom',...
                'YAxisLocation','right',...
                'XLim', [-6 10]...
           );
        ylabel('Normalized RMS');

        ax2 = axes('Position',get(ax1, 'Position'),...
                   'YAxisLocation','left',...
                   'Color','none',...
                   'XLim', [-6 10],...
                   'XColor','k','YColor','k');
        hold on
        % color coding for the posterior probabilities of the 3 locations
        stnc = 'r' ; prec = 'b'; postc = 'g' ;

        plot(ax2, depsort, pstn(depind), stnc, depsort, ppre(depind), prec,depsort, ppost(depind), postc, 'LineWidth', 2)
        set(ax1, 'XLim', [get(ax2, 'XLim')]);
        set(ax2, 'XDir', 'reverse', 'YLim', [0 1.5], 'YTick', [0 .5 1])
        lh = legend('In STN Probability', 'Before STN Probability', 'After STN Probability','Location', 'NorthWest');
        set(lh, 'Color', 'none');
        axes(ax1);
        arrpos = get(ax1, 'YLim');
  
        ar1 = line([stnin stnin], [0 arrpos(2)-0.5], 'Color', 'r', 'LineStyle', '--');
        ar2 = line([stnout stnout], [0 arrpos(2)-0.5],'Color', 'r', 'LineStyle', '--');
        axes(ax2);
        xlabel('EDT')
        ylabel('Probability')
    end
    