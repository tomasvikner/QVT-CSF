function MD = set_metadata(MD)
% Plot colors and LineWidth
c = linspecer(2);
MD.c1 = c(1, :); % plor color CSF 
MD.c2 = c(2, :); % plot color CBF 
MD.PLW = 2.4; % Waveform plot linewidth 
MD.MAXLAG = 250; % max CBF to CSF lag (ms) 
MD.wfState = 'CnoB'; % set initial CSF WF state 
MD.dotSize = 25;
MD.coregTtarget = 'cnob'; % or 'bcube';
MD.PATCH = true; % patch SD on waveform plot, otherwise multiple WFs 
end