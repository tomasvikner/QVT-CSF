function vMean = load_vMean(newDIM, MD, IDX)

% Just the blood flow data (TEMP: note scaling 1000)
disp('Computing time-averaged data')
vx = h5read(fullfile(MD.directory, MD.infile),'/bphx')/1000;
vy = h5read(fullfile(MD.directory, MD.infile),'/bphy')/1000;
vz = h5read(fullfile(MD.directory, MD.infile),'/bphz')/1000;
V(:,:,:,1) = single(mean(vx, 4)); % time-average here
V(:,:,:,2) = single(mean(vy, 4));
V(:,:,:,3) = single(mean(vz, 4));
vMean = zeros(newDIM(1),newDIM(2),newDIM(3),3,'single');
for n = 1:3
    vMean(:,:,:,n) = V(IDX.start(1):IDX.end(1),IDX.start(2):IDX.end(2),IDX.start(3):IDX.end(3),n);
end

end