function pindex = getBranchIndices(dcm_obj, branchList)
info = getCursorInfo(dcm_obj);
if isempty(info)
    pindex = 1;
    return;
end
ptList = reshape([info.Position], 3, []).';
pindex = zeros(size(ptList, 1), 1);
xyz = double(branchList(:, 1:min(3, size(branchList, 2))));
for n = 1:size(ptList, 1)
    d2 = sum((xyz - ptList(n, :)).^2, 2);
    [~, idxMin] = min(d2);
    pindex(n) = idxMin;
end
end
