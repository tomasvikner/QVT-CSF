function manualCoreg_Callback(~, ~, ~)
global dcm_obj CSFSEG CBFSEG branchList

disp('Called manual 2D point-point Co-Reg')
pindex = getBranchIndices(dcm_obj, branchList);
imdim = sqrt(numel(CBFSEG(1, :)));

moving = reshape(CSFSEG.cube(pindex, :), imdim, imdim);
fixed  = reshape(CSFSEG.mcsf(pindex, :), imdim, imdim); % Using mag only

figure(11); imshow(moving, []); [xm, ym] = ginput(1);
figure(12); imshow(fixed,  []); [xf, yf] = ginput(1);

dx = xf - xm; dy = yf - ym;
translated = imtranslate(moving, [dx, dy], 'OutputView', 'same');
figure(13); imshow(translated, []);

CSFSEG.cube(pindex, :)  = translated(:);
CSFSEG.bcube(pindex, :) = translated(:) > graythresh(translated(:));
CSFSEG.coregTrack(pindex) = 1;

set(dcm_obj, 'UpdateFcn', @myupdatefcn_all);
close([11 12 13]);
updateWaveforms('dcube');

function manualCoreg_Callback_3D(~, ~, ~)
global dcm_obj CSFSEG CBFSEG branchList

disp('Called manual 3D point-point Co-Reg')
pindex = getBranchIndices(dcm_obj, branchList);
imdim = sqrt(numel(CBFSEG(1,:)));

moving = reshape(CSFSEG.cube(pindex,:), imdim, imdim);
fixed  = reshape(CSFSEG.mcsf(pindex,:), imdim, imdim);

% Get multiple matching points
figure(11); imshow(moving, []); [xm, ym] = ginput(5);
figure(12); imshow(fixed,  []); [xf, yf] = ginput(5);

% Estimate and apply rigid 2D transformation
tform = fitgeotrans([xm ym], [xf yf], 'nonreflectivesimilarity');
translated = imwarp(moving, tform, 'OutputView', imref2d(size(moving)));

% Display result
figure(13); imshow(translated, []);

% Update segmentation
CSFSEG.cube(pindex,:)  = translated(:);
CSFSEG.bcube(pindex,:) = translated(:) > graythresh(translated(:));
CSFSEG.coregTrack(pindex) = 1;

set(dcm_obj, 'UpdateFcn', @myupdatefcn_all);
close([11 12 13]);
updateWaveforms('dcube');