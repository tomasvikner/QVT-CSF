
function [CBF, CSF, cbfwf, csfwf] = interpCoupling(METADATA, CBF, CSF, cbfwf, csfwf)
% Interpolate in time. xo length follows cbfwf columns (actual TR samples).
iframes = 1000;
if isfield(METADATA, 'iframes') && ~isempty(METADATA.iframes) && isnumeric(METADATA.iframes)
    iframes = max(2, round(double(METADATA.iframes(1))));
end

% At least one column (time sample) and one row (point) so interp1 sees a non-empty matrix
if size(cbfwf, 2) < 1
    cbfwf = zeros(max(1, size(cbfwf, 1)), 1);
end
if size(cbfwf, 1) < 1
    cbfwf = zeros(1, size(cbfwf, 2));
end
nData = size(cbfwf, 2);

if size(csfwf, 2) < nData
    csfwf(:, end+1:nData) = 0;
elseif size(csfwf, 2) > nData
    csfwf = csfwf(:, 1:nData);
end
if size(csfwf, 1) < 1
    csfwf = zeros(1, nData);
end

CBF = CBF(:)';
CSF = CSF(:)';
if numel(CBF) < nData
    CBF(1, nData) = 0;
elseif numel(CBF) > nData
    CBF = CBF(1:nData);
end
if numel(CSF) < nData
    CSF(1, nData) = 0;
elseif numel(CSF) > nData
    CSF = CSF(1:nData);
end

if size(cbfwf, 2) < nData
    cbfwf(:, end+1:nData) = 0;
elseif size(cbfwf, 2) > nData
    cbfwf = cbfwf(:, 1:nData);
end

if nData < 2
    cbfwf = [cbfwf(:, 1) cbfwf(:, 1)];
    csfwf = [csfwf(:, 1) csfwf(:, 1)];
    if numel(CBF) >= 1
        CBF = [CBF(1) CBF(1)];
    else
        CBF = [0 0];
    end
    if numel(CSF) >= 1
        CSF = [CSF(1) CSF(1)];
    else
        CSF = [0 0];
    end
    nData = 2;
end

xo = linspace(1, nData, nData);
xq = linspace(1, nData, iframes);
CBF = interp1(xo, CBF, xq, 'pchip');
CSF = interp1(xo, CSF, xq, 'pchip');
cbfwf = interp1(xo, cbfwf', xq, 'pchip');
csfwf = interp1(xo, csfwf', xq, 'pchip');
end
