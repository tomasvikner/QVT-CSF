function IDX = calc_crop(MAG)

% Auto crop images to save memory 
SUMnumA = squeeze(sum(sum(MAG,1),2)); %1D axial projection
SUMnumS = squeeze(sum(sum(MAG,1),3))'; %1D sagittal projection
SUMnumC = squeeze(sum(sum(MAG,2),3)); %1D coronal projection

% Chop off edges of projections (usually noisy data)
SUMnumA(1:3) = 0; SUMnumA(end-2:end) = 0;
SUMnumS(1:3) = 0; SUMnumS(end-2:end) = 0;
SUMnumC(1:3) = 0; SUMnumC(end-2:end) = 0;

% Normalize values (from 0-1)
SUMnumC = rescale(SUMnumC,'InputMin',min(SUMnumC),'InputMax',max(SUMnumC)); 
BIN = SUMnumC>0.25; % Find where projection crosses thresh value of 0.25 
[~,IDX.start(1)] = max(BIN,[],1); %get first thresh crossing
[~,IDX.end(1)] = max(flipud(BIN),[],1); 
IDX.end(1) = matrix(1) - IDX.end(1) + 1; %get last thresh crossing

SUMnumS = rescale(SUMnumS,'InputMin',min(SUMnumS),'InputMax',max(SUMnumS)); 
BIN = SUMnumS>0.25; % Find where projection crosses thresh value of 0.25 
[~,IDX.start(2)] = max(BIN,[],1); %get first thresh crossing
[~,IDX.end(2)] = max(flipud(BIN),[],1); 
IDX.end(2) = matrix(2) - IDX.end(2) + 1; %get last thresh crossing

SUMnumA = rescale(SUMnumA,'InputMin',min(SUMnumA),'InputMax',max(SUMnumA)); 
BIN = SUMnumA>0.25; % Find where projection crosses thresh value of 0.25 
[~,IDX.start(3)] = max(BIN,[],1); %get first thresh crossing
[~,IDX.end(3)] = max(flipud(BIN),[],1); 
IDX.end(3) = matrix(3) - IDX.end(3) + 1; %get last thresh crossing

end