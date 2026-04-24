function cyl = center_cylinder_mask(sz, radiusVox)
%CENTER_CYLINDER_MASK  Solid z-axis cylinder through FOV xy-center (no intensity threshold).
%   sz = [nx ny nz]; radiusVox is Euclidean radius in the x–y plane (voxels, default 20).
%   Every slice gets the same filled disk (inclusive boundary).

if nargin < 2
    radiusVox = [];
end
if isempty(radiusVox)
    radiusVox = 20;
end
sz = double(sz(:).');
if numel(sz) < 3
    cyl = false;
    return;
end
nx = sz(1);
ny = sz(2);
nz = sz(3);
cx = (nx + 1) / 2;
cy = (ny + 1) / 2;
[xg, yg] = ndgrid(1:nx, 1:ny);
rsq = double(radiusVox) ^ 2;
disk = ((xg - cx) .^ 2 + (yg - cy) .^ 2) <= rsq;
cyl = false(nx, ny, nz);
for z = 1:nz
    cyl(:, :, z) = disk;
end

end
