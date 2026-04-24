function segmentOut = preview_centerline_segmentation(cdMasked, segmentFallback, MAG)
%PREVIEW_CENTERLINE_SEGMENTATION  Blocking QC: threshold + cylinder CC, then continue pipeline.
%   Blocks loadHDF5 until you click "Export and continue" or close the window.
%   Close without Export → uses segmentFallback (automatic segment_angiogram result).
%
%   segmentOut = preview_centerline_segmentation(cdMasked, segmentFallback, MAG)
%
%   cdMasked: masked CD magnitude. segmentFallback: initial auto mask (slider hint + default).
%   For display, geometric cylinder is OR-ed in (solid column); segmentOut never includes
%   that cylinder — only CC-filtered threshold mask (same as Export body).

if nargin < 3
    MAG = [];
end

tag = 'QVT_CenterlineSegPreview';
old = findall(0, 'Type', 'figure', 'Tag', tag);
if ~isempty(old)
    delete(old);
end

vol = single(cdMasked);
[s1, s2, s3] = size(vol);
z0 = max(1, round(s3 / 2));
magFrac = 0.07;
cylR = 20;
cylConn = 26;

if isempty(MAG) || ~isequal(size(MAG), size(vol))
    tissueMask = true(s1, s2, s3);
else
    mag = single(MAG);
    mx = max(mag(:));
    if mx <= 0
        tissueMask = true(s1, s2, s3);
    else
        tissueMask = mag > (magFrac * mx);
        tissueMask = imfill(tissueMask, 'holes');
    end
end

maxV = double(max(vol(:)));
if maxV <= 0
    maxV = 1;
end

vals = vol(segmentFallback(:) & tissueMask(:));
if isempty(vals)
    tRel0 = 0.2;
else
    tRel0 = median(double(vals)) / maxV;
end
tRel0 = min(max(tRel0, 0.02), 0.98);

f = figure('IntegerHandle', 'off', ...
    'HandleVisibility', 'callback', ...
    'NumberTitle', 'off', ...
    'WindowStyle', 'modal', ...
    'Name', 'QVT: adjust threshold — Export and continue (or close for auto mask)', ...
    'Tag', tag, ...
    'Position', [80 80 920 700]);

t = tiledlayout(f, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact', ...
    'Units', 'normalized', 'Position', [0.04 0.20 0.92 0.76]);
axAxOv = nexttile(t, 1);
axMipOv = nexttile(t, 2);
axAxBin = nexttile(t, 3);
axMipBin = nexttile(t, 4);

thrSlider = uicontrol('Parent', f, 'Style', 'slider', 'Units', 'normalized', ...
    'Position', [0.12 0.13 0.76 0.035], 'Min', 0.001, 'Max', 0.999, 'Value', tRel0, ...
    'TooltipString', 'Global threshold: fraction of max angiogram (inside tissue mask)');
thrLabel = uicontrol('Parent', f, 'Style', 'text', 'Units', 'normalized', ...
    'Position', [0.04 0.168 0.70 0.028], 'String', '', ...
    'BackgroundColor', get(f, 'Color'), 'HorizontalAlignment', 'left', ...
    'FontSize', 9);

zSlider = uicontrol('Parent', f, 'Style', 'slider', 'Units', 'normalized', ...
    'Position', [0.12 0.085 0.76 0.035], 'Min', 1, 'Max', s3, 'Value', z0, ...
    'SliderStep', [1 / max(s3 - 1, 1), 10 / max(s3 - 1, 1)], ...
    'TooltipString', 'Axial slice index');
uicontrol('Parent', f, 'Style', 'text', 'Units', 'normalized', ...
    'Position', [0.12 0.048 0.76 0.028], 'String', 'Axial slice (z)', ...
    'BackgroundColor', get(f, 'Color'), 'HorizontalAlignment', 'center');

uicontrol('Parent', f, 'Style', 'pushbutton', 'Units', 'normalized', ...
    'Position', [0.76 0.168 0.20 0.034], 'String', 'Export and continue', ...
    'TooltipString', 'Apply current CC-filtered mask and resume loading (also saves QVT_segPreview)', ...
    'Callback', @onExport);

uicontrol('Parent', f, 'Style', 'text', 'Units', 'normalized', ...
    'Position', [0.04 0.005 0.92 0.038], 'String', ...
    'Pipeline waits here. Close window to keep the automatic mask instead.', ...
    'BackgroundColor', get(f, 'Color'), 'HorizontalAlignment', 'left', ...
    'FontSize', 9);

set(f, 'UserData', struct('result', [], 'fallback', segmentFallback, 'exported', false));
set(f, 'CloseRequestFcn', @onCloseFigure);

    function seg = segFromThr(tr)
        seg = (vol >= single(tr * maxV)) & tissueMask;
    end

    function seg = segThrCylinder(tr)
        seg = filter_segment_touching_cylinder(segFromThr(tr), cylR, cylConn);
    end

    function drawAll(z, tr)
        z = max(1, min(s3, round(z)));
        seg = segThrCylinder(tr);
        cylG = center_cylinder_mask([s1, s2, s3], cylR);
        segVis = seg | cylG;
        cdN = mat2gray(double(vol));
        A = cdN(:, :, z);
        M = segVis(:, :, z);
        R = min(A + 0.5 * double(M), 1);
        G = A .* (1 - 0.35 * double(M));
        B = A .* (1 - 0.35 * double(M));
        imagesc(axAxOv, cat(3, R, G, B));
        axis(axAxOv, 'image');
        axis(axAxOv, 'off');
        title(axAxOv, sprintf('Axial z=%d/%d — angio + mask', z, s3), 'FontSize', 10);

        imagesc(axAxBin, M);
        axis(axAxBin, 'image');
        axis(axAxBin, 'off');
        colormap(axAxBin, [0 0 0; 1 0.35 0.1]);
        title(axAxBin, 'Axial: CC mask | solid cyl (geom)', 'FontSize', 10);

        mipA = squeeze(max(vol, [], 3));
        mipS = squeeze(max(single(segVis), [], 3));
        g = mat2gray(double(mipA));
        R2 = min(g + 0.55 * double(mipS), 1);
        G2 = g .* (1 - 0.4 * double(mipS));
        B2 = g .* (1 - 0.4 * double(mipS));
        imagesc(axMipOv, cat(3, R2, G2, B2));
        axis(axMipOv, 'image');
        axis(axMipOv, 'off');
        title(axMipOv, 'MIP z (max along dim 3): angio + mask', 'FontSize', 10);

        imagesc(axMipBin, mipS);
        axis(axMipBin, 'image');
        axis(axMipBin, 'off');
        colormap(axMipBin, [0 0 0; 1 0.35 0.1]);
        title(axMipBin, 'MIP z: CC mask | solid cyl (geom)', 'FontSize', 10);

        trv = tr * maxV;
        set(thrLabel, 'String', sprintf( ...
            'Threshold = %.4f × max (masked) = %.5g  |  CC mask voxels: %d', ...
            tr, trv, nnz(seg)));
        drawnow limitrate;
    end

    function onThr(~, ~)
        tr = get(thrSlider, 'Value');
        z = round(get(zSlider, 'Value'));
        drawAll(z, tr);
    end

    function onZ(~, ~)
        tr = get(thrSlider, 'Value');
        z = round(get(zSlider, 'Value'));
        drawAll(z, tr);
    end

    function onExport(~, ~)
        tr = get(thrSlider, 'Value');
        mask = segThrCylinder(tr);
        ud = get(f, 'UserData');
        ud.result = mask;
        ud.exported = true;
        set(f, 'UserData', ud);
        assignin('base', 'QVT_segPreview', mask);
        disp('QVT: continuing pipeline with exported CC-filtered mask (QVT_segPreview).');
        uiresume(f);
    end

    function onCloseFigure(src, ~)
        ud = get(src, 'UserData');
        if isempty(ud.result)
            ud.result = ud.fallback;
            if ~ud.exported
                disp('QVT: continuing pipeline with automatic segmentation (window closed, no Export).');
            end
        end
        set(src, 'UserData', ud);
        uiresume(src);
    end

set(thrSlider, 'Callback', @onThr);
set(zSlider, 'Callback', @onZ);

drawAll(z0, tRel0);

uiwait(f);

if isgraphics(f)
    ud = get(f, 'UserData');
    segmentOut = ud.result;
    if isempty(segmentOut)
        segmentOut = ud.fallback;
    end
    set(f, 'CloseRequestFcn', []);
    delete(f);
else
    segmentOut = segmentFallback;
end

segmentOut = logical(segmentOut);

end
