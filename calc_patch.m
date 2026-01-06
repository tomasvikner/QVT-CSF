function [X, Y] = calc_patch(WF, cardiacCycle)   

    sdwf = std(WF, [], 2)';
    upper_bound = WF + sdwf;
    lower_bound = WF - sdwf;
    X = [cardiacCycle(:)', fliplr(cardiacCycle(:)')];
    Y = [upper_bound(:)', fliplr(lower_bound(:)')];

end