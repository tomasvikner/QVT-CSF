
function [CBF, CSF, cbfwf, csfwf] = interpCoupling(METADATA, CBF, CSF, cbfwf, csfwf)
nframes = METADATA.nframes; % Original frames
iframes = METADATA.iframes; % Desired output resolution
xo = linspace(1, nframes, nframes);
xq = linspace(1, nframes, iframes);
CBF = interp1(xo, CBF', xq, 'pchip'); % mean point
CSF = interp1(xo, CSF', xq, 'pchip');
cbfwf = interp1(xo, cbfwf', xq, 'pchip');
csfwf = interp1(xo, csfwf', xq, 'pchip'); % all points
end
