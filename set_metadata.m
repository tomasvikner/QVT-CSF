function MD = set_metadata(MD)
% Plot colors and LineWidth
c = linspecer(2);
MD.c1 = c(1, :); % plor color CSF 
MD.c2 = c(2, :); % plot color CBF 
MD.PLW = 2.4; % Waveform plot linewidth 
MD.MAXLAG = 250; % max CBF to CSF lag (ms) 
MD.wfState = 'CnoB'; % set initial CSF WF state 
MD.dotSize = 25;
MD.coregTarget = 'cnob'; % or 'bcube'
MD.PATCH = true; % patch SD on waveform plot, otherwise multiple WFs 
% High-res cardiac phase grid for interpCoupling / plots (standard = 1000 samples)
MD.iframes = 1000;
if ~isfield(MD, 'RECONRES') || isempty(MD.RECONRES)
    MD.RECONRES = 0.75; % mm isotropic default only if not set (e.g. after loadHDF5 resampling)
end
if ~isfield(MD, 'nframes') || isempty(MD.nframes)
    MD.nframes = 20;
end
if ~isfield(MD, 'INTERP') || isempty(MD.INTERP)
    MD.INTERP = false; % true: HDF5 load uses isotropic resampling in loadHDF5
end
end