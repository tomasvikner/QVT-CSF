function pindex = getBranchIndices(dcm_obj, branchList)
info = getCursorInfo(dcm_obj);
ptList = reshape([info.Position], 3, []).';
pindex = zeros(size(ptList, 1), 1);

for n = 1:size(ptList,1)
    idx = all(bsxfun(@eq, branchList(:,1:3), ptList(n,:)), 2);
    pindex(n) = find(idx, 1);
end
