% --- help function for updating waveforms from sliders/mansegs
function [pindex, index_range] = get_pindex(dcm_obj, branchList)
info_struct = getCursorInfo(dcm_obj);
ptList = [info_struct.Position];
ptList = reshape(ptList,[3,numel(ptList)/3])';
pindex = zeros(size(ptList,1),1);
for n = 1:size(ptList,1)
    xIdx = find(branchList(:,1) == ptList(n,1));
    yIdx = find(branchList(xIdx,2) == ptList(n,2));
    zIdx = find(branchList(xIdx(yIdx),3) == ptList(n,3));
    pindex(n) = xIdx(yIdx(zIdx)); %#ok<*FNDSB>
end
bnum = branchList(pindex,4);
Logical_branch = branchList(:,4) ~= bnum;
index_range = pindex;
index_range(index_range<1) = [];
index_range(index_range>size(branchList,1)) = [];
index_range(Logical_branch(index_range)) = [];
end