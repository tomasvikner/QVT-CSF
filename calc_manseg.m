function mask = calc_manseg(dcm_obj, branchList, CBFSEG)

disp('Called manual 2D segmentation')
pindex = getBranchIndices(dcm_obj, branchList);
imdim = sqrt(numel(CBFSEG(1,:)));

% moving = reshape(CSFSEG.cube(pindex,:), imdim, imdim);
moving = reshape(CSFSEG.scsf(pindex,:), imdim, imdim); % Manual seg on cnob, not cube 

figure(11); imshow(moving, []);
disp('Draw region with freehand tool. Double-click to finish.');
clim([min(moving(:)), 0.5 * max(moving(:))]);

h = imfreehand;
mask = createMask(h);
close(11);

end