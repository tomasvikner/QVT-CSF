function [amp, cdv] = calc_metrics(CSF, CBF)
amp = [];
amp.CSF = max(CSF) - min(CSF);
amp.CBF = max(CBF) - min(CBF);
cdv = [];
ccsf = cumtrapz(CSF - mean(CSF)) * 1e-3; % with current interp? 
ccbf = cumtrapz(CBF - mean(CBF)) * 1e-3;
cdv.CSF = max(ccsf) - min(ccsf);
cdv.CBF = max(ccbf) - min(ccbf);
end
