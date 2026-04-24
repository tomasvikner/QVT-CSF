function segment = filter_segment_touching_cylinder(segment, radiusVox, conn3d)
%FILTER_SEGMENT_TOUCHING_CYLINDER  Keep only 3-D CCs that intersect a central cylinder.
%   Cylinder axis = z (column through geometric center); solid filled disk in x–y
%   (see center_cylinder_mask — geometric only, no magnitude threshold).
%
%   segment = filter_segment_touching_cylinder(segment)
%   segment = filter_segment_touching_cylinder(segment, radiusVox, conn3d)
%     radiusVox — default 20
%     conn3d    — connectivity for bwconncomp (default 26)

if nargin < 2 || isempty(radiusVox)
    radiusVox = 20;
end
if nargin < 3 || isempty(conn3d)
    conn3d = 26;
end

segment = logical(segment);
sz = size(segment);
if numel(sz) < 3
    return;
end

cyl = center_cylinder_mask(sz, radiusVox);

cc = bwconncomp(segment, conn3d);
if cc.NumObjects == 0
    segment = false(sz);
    return;
end

segOut = false(sz);
for k = 1:cc.NumObjects
    idx = cc.PixelIdxList{k};
    if any(cyl(idx))
        segOut(idx) = true;
    end
end
segment = segOut;

end
