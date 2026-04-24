% --- help function for updating waveforms from sliders/mansegs
function [pindex, index_range, pinds] = get_pindex(dcm_obj, branchList)
info_struct = getCursorInfo(dcm_obj);
if isempty(info_struct)
    pindex = 1;
    index_range = 1;
    pinds = 1;
    return;
end
ptList = [info_struct.Position];
ptList = reshape(ptList, [3, numel(ptList) / 3])';

xyz = double(branchList(:, 1:min(3, size(branchList, 2))));
pindex = zeros(size(ptList, 1), 1);
for n = 1:size(ptList, 1)
    d2 = sum((xyz - ptList(n, :)).^2, 2);
    [~, idxMin] = min(d2);
    pindex(n) = idxMin;
end

bnum = branchList(pindex(1), 4);
Logical_branch = branchList(:, 4) ~= bnum;
index_range = pindex(:);
index_range(index_range < 1) = [];
index_range(index_range > size(branchList, 1)) = [];
index_range(Logical_branch(index_range)) = [];
if isempty(index_range)
    index_range = pindex(1);
end

pside = index_range;
pside(pside == pindex(1)) = [];
if size(pside, 2) > size(pside, 1)
    pside = pside';
end
pinds = [pindex(1); pside];
end
