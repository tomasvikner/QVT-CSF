function [rmax, mlag] = waveformCoupling(wf1, wf2, maxlag)

    % Wrap the waveform 
    wf1 = [wf1(:, end-maxlag+1:end), wf1, wf1(:, 1:maxlag)];
    wf2 = [wf2(:, end-maxlag+1:end), wf2, wf2(:, 1:maxlag)];
    wf1 = zscore(wf1);
    wf2 = - zscore(wf2);
    
    % Cross-correlation
    [rvals, lags] = xcorr(wf2, wf1, maxlag, 'coeff'); 
    [~, mi] = max(abs(rvals));
    mlag = lags(mi);
    rmax = - rvals(mi);

end