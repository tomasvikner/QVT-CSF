function [CBF, CSF, cbfwf, csfwf] = calc_waveforms(CLVALS, flowCSF, pinds, wfState)

cbfwf = CLVALS.flowTR(pinds, :);
if strcmp(wfState, 'PC-1')
    csfwf = flowCSF.PC1.mean(pinds, :); 
elseif strcmp(wfState, 'CnoB')
    csfwf = flowCSF.cnob.mean(pinds, :); 
elseif strcmp(wfState, 'CUBE')
    csfwf = flowCSF.bcube.mean(pinds, :); 
end

zeroRows = all(csfwf == 0, 2); 
nanRows = all(isnan(csfwf), 2);
emptyRows = find(zeroRows | nanRows);
csfwf(emptyRows, :) = []; 
cbfwf(emptyRows, :) = []; 

ncsf = size(csfwf, 1);
nT = size(cbfwf, 2);
if nT < 1
    nT = 1;
end
if ncsf == 0
    CSF = zeros(1, nT);
elseif ncsf == 1
    CSF = csfwf;
else
    CSF = median(csfwf);
end
CBF = median(cbfwf);

if isempty(CSF)
    disp('CSF waveform empty!')
    disp(size(cbfwf))
end

if isempty(CBF)
    disp('CBF waveform empty!')
end