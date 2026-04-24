function [StructVols, MAG, vMean, MD, vel4dResampled] = resample_pcvipr_to_isotropic(StructVols, MAG, vel4dNative, MD)
%RESAMPLE_PCVIPR_TO_ISOTROPIC  After crop: isotropic 0.5 mm grid + nFrames temporal resampling.
%   FoV default 4 x 24 x 24 cm (40 x 240 x 240 mm) aligned with array dims 1:3 of MAG.
%   Output grid: round(fov_mm / 0.5) => [80 480 480] for that FoV.
%   Temporal: pick MD.nframes_target (default 20) cardiac phases by nearest native index (no interp1).
%   Spatial: imresize3 only (per frame, linear).
%
%   Updates MD only with scalar/grid metadata: RECONRES, nframes, nativeSpacing_mm, fov_mm,
%   targetIsotropic_mm, reconSize_nn, nframes_native. Resampled 4-D CBF+CSF velocities are
%   returned in vel4dResampled (vxf,vyf,vzf,cxf,cyf,czf) — not attached to MD.

target_mm = 0.5;
if isfield(MD, 'targetIsotropic_mm') && ~isempty(MD.targetIsotropic_mm)
    target_mm = double(MD.targetIsotropic_mm(1));
end

nFramesOut = 20;
if isfield(MD, 'nframes_target') && ~isempty(MD.nframes_target)
    nFramesOut = round(double(MD.nframes_target(1)));
end
nFramesOut = max(2, nFramesOut);

fov_mm = [40, 240, 240];
if isfield(MD, 'fov_mm') && numel(MD.fov_mm) >= 3
    fov_mm = double(MD.fov_mm(1:3));
end

sz = size(MAG);
if numel(sz) < 3
    error('resample_pcvipr_to_isotropic:MAG', 'MAG must be 3-D after crop.');
end

native_mm = fov_mm(:) ./ sz(:);

if isfield(MD, 'interpGrid_nn') && ~isempty(MD.interpGrid_nn) && numel(MD.interpGrid_nn) >= 3
    outSz = round(double(MD.interpGrid_nn(1:3)));
else
    outSz = round(fov_mm(:) ./ target_mm);
end
outSz = max(outSz, 2);

MD.nframes_native = size(vel4dNative.bphx, 4);
MD.nativeSpacing_mm = native_mm.';
MD.fov_mm = fov_mm(:).';
MD.targetIsotropic_mm = target_mm;
MD.reconSize_nn = outSz(:).';
MD.RECONRES = target_mm;
MD.nframes = nFramesOut;

% --- 3-D magnitude / std volumes ---
fns = fieldnames(StructVols);
for k = 1:numel(fns)
    fn = fns{k};
    V = StructVols.(fn);
    if ~isequal(size(V), sz)
        continue;
    end
    StructVols.(fn) = imresize3(single(V), outSz, 'linear');
end
MAG = StructVols.mcbf;

% --- 4-D velocity: CBF (blood) then CSF — temporal then spatial ---
nmB = {'bphx', 'bphy', 'bphz'};
nmC = {'cphx', 'cphy', 'cphz'};
vel4dResampled = struct();
for ii = 1:3
    V = single(vel4dNative.(nmB{ii}));
    V = resample_time_4d(V, nFramesOut);
    V = resample_space_4d(V, outSz);
    switch ii
        case 1, vel4dResampled.vxf = V;
        case 2, vel4dResampled.vyf = V;
        case 3, vel4dResampled.vzf = V;
    end
end
for ii = 1:3
    V = single(vel4dNative.(nmC{ii}));
    V = resample_time_4d(V, nFramesOut);
    V = resample_space_4d(V, outSz);
    switch ii
        case 1, vel4dResampled.cxf = V;
        case 2, vel4dResampled.cyf = V;
        case 3, vel4dResampled.czf = V;
    end
end

% Time-averaged CBF (blood) velocity — same ÷1000 convention as load_vMean
vMean = zeros(outSz(1), outSz(2), outSz(3), 3, 'single');
vMean(:, :, :, 1) = mean(vel4dResampled.vxf, 4, 'omitnan');
vMean(:, :, :, 2) = mean(vel4dResampled.vyf, 4, 'omitnan');
vMean(:, :, :, 3) = mean(vel4dResampled.vzf, 4, 'omitnan');

end

function V = resample_time_4d(V, nOut)
% Nearest-phase indexing only (avoids interp1 over nx*ny*nz rows — huge memory).
[nx, ny, nz, nT] = size(V);
if nT == nOut
    return;
end
if nT < 2
    V = repmat(V, [1 1 1 nOut]);
    return;
end
tq = linspace(1, nT, nOut);
idx = min(nT, max(1, round(tq)));
V = V(:, :, :, idx);
V = single(V);
end

function V = resample_space_4d(V, outSz)
[nx, ny, nz, nT] = size(V);
Vo = zeros([outSz(:)', nT], 'like', V);
for t = 1:nT
    Vo(:, :, :, t) = imresize3(V(:, :, :, t), outSz, 'linear');
end
V = single(Vo);
end
