function [segment, cdMasked] = segment_angiogram_centerline(cdImg, MAG, varargin)
%SEGMENT_ANGIOGRAM_CENTERLINE  Global vessel mask for centerline extraction.
%   cdImg: 3-D complex-difference–style magnitude (see calc_cd_timeavg), not a temporal MIP.
%   Constrains by magnitude-based tissue mask, slidingThreshold, then cleanup.
%
%   [segment, cdMasked] = segment_angiogram_centerline(cdImg, MAG)
%   Optional name-value pairs:
%     'MagFrac'   — lower MAG threshold as fraction of max(MAG) (default 0.07)
%     'AreaFrac'  — bwareaopen minimum area = max(round(frac*foreground), AbsMin)
%     'AreaMin'   — absolute minimum component size in voxels (default 150)
%     'Conn'      — connectivity for bwareaopen (default 6)
%     'CylinderRadius' — voxels; z-axis cylinder at FOV xy-center (default 20)
%     'CylinderConn'   — bwconncomp connectivity for cylinder filter (default 26)

p = inputParser;
addParameter(p, 'MagFrac', 0.07, @(x) isnumeric(x) && isscalar(x) && x > 0 && x < 1);
addParameter(p, 'AreaFrac', 0.012, @(x) isnumeric(x) && isscalar(x) && x > 0 && x < 1);
addParameter(p, 'AreaMin', 150, @(x) isnumeric(x) && isscalar(x) && x >= 0);
addParameter(p, 'Conn', 6, @(x) isnumeric(x) && ismember(x, [6 18 26]));
addParameter(p, 'CylinderRadius', 20, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'CylinderConn', 26, @(x) isnumeric(x) && ismember(x, [6 18 26]));
parse(p, varargin{:});
magFrac = p.Results.MagFrac;
areaFrac = p.Results.AreaFrac;
areaMin = p.Results.AreaMin;
conn = p.Results.Conn;
cylR = p.Results.CylinderRadius;
cylConn = p.Results.CylinderConn;

cdImg = single(cdImg);
MAG = single(MAG);

magMax = max(MAG(:));
if magMax <= 0
    tissueMask = true(size(MAG));
else
    tissueMask = MAG > (magFrac * magMax);
    tissueMask = imfill(tissueMask, 'holes');
end

cdMasked = cdImg .* tissueMask;

step = 0.001;
UPthresh = 0.8;
SMf = 10;
shiftHM_flag = 1;
medFilt_flag = 1;
[~, segment] = slidingThreshold(cdMasked, step, UPthresh, SMf, shiftHM_flag, medFilt_flag);

segment = logical(segment) & tissueMask;
% Keep only thresholded components that meet a central z-axis cylinder (FOV center)
segment = filter_segment_touching_cylinder(segment, cylR, cylConn);
segment = imclearborder(segment);

S = sum(segment(:));
areaThresh = max(round(S * areaFrac), areaMin);
segment = bwareaopen(segment, areaThresh, conn);

CC = bwconncomp(segment, conn);
if CC.NumObjects == 0
    return;
end
if CC.NumObjects > 1
    np = cellfun(@numel, CC.PixelIdxList);
    [sortedA, ord] = sort(np, 'descend');
    cumA = cumsum(sortedA);
    keep = cumA <= 0.98 * sum(np);
    if ~any(keep)
        keep(1) = true;
    end
    keepIdx = ord(keep);
    if isempty(keepIdx)
        keepIdx = ord(1);
    end
    segNew = false(size(segment));
    for k = 1:numel(keepIdx)
        segNew(CC.PixelIdxList{keepIdx(k)}) = true;
    end
    segment = segNew;
end

end
