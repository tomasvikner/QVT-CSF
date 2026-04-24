function diam_val = calc_diam(segment)
    % New with ratios of areas. Ratio of smallest inner circle over
    % largest encompassing outer circle (assume circular area). Measure of
    % circularity of vessel (ratio =1 is circle,ratio<1 is irregular shape)
    D = bwdist(~segment); %euclidean distance transform
    Rin = max(D(:)); %distance from center to closest non-zero entry
    [xLoc,yLoc] = find(bwperim(segment)); %get perimeter
    D = pdist2([xLoc,yLoc],[xLoc,yLoc]); %distance b/w perimeter points
    Rout = max(D(:))/2; %radius of largest outer circle
    diam_val = Rin^2/Rout^2; %ratio of areas
    diam_val(diam_val==inf) = 0;
end