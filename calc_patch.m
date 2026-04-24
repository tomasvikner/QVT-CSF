function [X, Y] = calc_patch(WF, cardiacCycle)
%CALC_PATCH  Envelope (mean ± std across traces) for PATCH waveform shading.
%   WF is nTime x nTraces after interpCoupling; cardiacCycle has length nTime.

nt = numel(cardiacCycle);
WF = double(WF);
if isempty(WF)
    X = [cardiacCycle(:)', fliplr(cardiacCycle(:)')];
    Y = zeros(size(X));
    return;
end

% Prefer time along rows
if size(WF, 1) ~= nt && size(WF, 2) == nt
    WF = WF.';
end
if size(WF, 1) ~= nt
    if size(WF, 1) > 0 && nt > 1
        tOld = linspace(1, size(WF, 1), size(WF, 1));
        tNew = linspace(1, size(WF, 1), nt);
        WF = interp1(tOld, WF, tNew, 'linear', 'extrap');
    else
        WF = repmat(WF(1, :), nt, 1);
    end
end

if size(WF, 2) <= 1
    m = WF(:, 1);
    sd = zeros(size(m));
else
    m = mean(WF, 2, 'omitnan');
    sd = std(WF, 0, 2, 'omitnan');
end
sd(~isfinite(sd)) = 0;
m(~isfinite(m)) = 0;

m = m(:)';
sd = sd(:)';
upper = min(m + sd, max(m) + eps);
lower = max(m - sd, min(m) - eps);
cc = cardiacCycle(:)';
X = [cc, fliplr(cc)];
Y = [upper, fliplr(lower)];

end
