function [] = save_point(info_struct, branchList, CLVALS, VcrossTR, CcrossTR, CC, CSFROI, flowCSF, WFPS, PointLabel, directory)

ptList = [info_struct.Position];
ptList = reshape(ptList,[3,numel(ptList)/3])';
pindex = zeros(size(ptList,1),1);

for n = 1:size(ptList,1)
    xIdx = find(branchList(:,1) == ptList(n,1));
    yIdx = find(branchList(xIdx,2) == ptList(n,2));
    zIdx = find(branchList(xIdx(yIdx),3) == ptList(n,3));
    pindex(n) = xIdx(yIdx(zIdx));
end

% Gives associated branch number if full branch point is wanted
bnum = branchList(pindex,4);
Logical_branch = branchList(:,4) ~= bnum;

% OUTPUT +/- points use this
index_range = pindex-2:pindex+2;

%removes outliers and points from other branches
index_range(index_range<1) = [];
index_range(index_range>size(branchList,1)) = [];
index_range(Logical_branch(index_range)) = [];

% disp('flowPulsatile (paramMap.m)')
% Time-resolved flow
flowPulsatile = flowPulsatile_val(index_range,:);
flowPulsatile = [flowPulsatile;mean(flowPulsatile,1);std(flowPulsatile,1)];

savePoint = [];

savePoint.TRCBF = VcrossTR;
savePoint.TRCSF = CcrossTR;
% savePoint.BMASK = Maskcross;
savePoint.CMASK = CC.madj;

savePoint.bnum = bnum;
savePoint.Logical_branch = Logical_branch;
savePoint.index_range = index_range;
savePoint.ptList = ptList;

savePoint.CBF = flowPulsatile;
savePoint.CSF.bics = flowCSF.bics.mean(index_range, :);
savePoint.CSF.bcsf = flowCSF.bcsf.mean(index_range, :);
savePoint.CSF.cnob = flowCSF.cnob.mean(index_range, :);
savePoint.CSF.full = flowCSF.full.mean(index_range, :);
savePoint.CSF.pc1 = flowCSF.PC1.mean(index_range, :);
DOMADJ = false;
if DOMADJ
    if sum(flowCSF.madj.mean(index_range, :)) > 0 %#ok<UNRCH>
        savePoint.CSF.madj = flowCSF.madj.mean(index_range, :);
    end
    savePoint.CSF.madj = flowCSF.madj.mean(index_range, :);
end

savePoint.CSFROI.bics = CSFROI.bics(index_range, :, :, :);
savePoint.CSFROI.bcsf = CSFROI.bcsf(index_range, :, :, :);
savePoint.CSFROI.cnob = CSFROI.cnob(index_range, :, :, :);
savePoint.CSFROI.flow = CSFROI.flow(index_range, :, :, :);

% Name of point label e.g. Left MCA as in scroll menu in GUI
savePoint.pointName = PointLabel;

% note - Waveform point save (struct) - THIS contains most important information to export!
savePoint.WFPS = WFPS;

% save([SavePath savePoint.pointName '.mat'], 'savePoint');
SavePath = [directory filesep 'PointsV4/'];
if ~exist(SavePath, 'dir')
    mkdir(SavePath)
end
save([SavePath savePoint.pointName '.mat'], 'savePoint');
disp(['Saved point @: ' [SavePath savePoint.pointName '.mat'] ]);
