function [VcrossTR, CcrossTR] = get_vCS_TR(Vplanes, pindex, imdim, METADATA)
%GET_VCS_TR  Time-resolved CBF/CSF cross-section velocities for one centerline point.
%   imdim     optional output side length (default = native plane side from Vplanes)
%   METADATA  optional; when set, temporal length uses METADATA.nframes (clamped to Vplanes)

    if nargin < 3 || isempty(imdim)
        imdim = [];
    end
    if nargin < 4
        METADATA = [];
    end

    [cbfStack, nFromCbf] = plane_stack_from_vplanes(Vplanes.CBF, pindex);
    [csfStack, nFromCsf] = plane_stack_from_vplanes(Vplanes.CSF, pindex);
    nfData = min(nFromCbf, nFromCsf);

    if ~isempty(METADATA) && isfield(METADATA, 'nframes') && ~isempty(METADATA.nframes)
        nframes = max(1, min(double(METADATA.nframes), nfData));
    else
        nframes = nfData;
    end

    cbfStack = cbfStack(:, 1:nframes);
    csfStack = csfStack(:, 1:nframes);

    VcrossTR = reshape_stack(cbfStack, 'CBF');
    CcrossTR = reshape_stack(csfStack, 'CSF');

    if isempty(imdim)
        imdim = size(VcrossTR, 1);
    end
    VcrossTR = resize_plane_stack(VcrossTR, imdim);
    CcrossTR = resize_plane_stack(CcrossTR, imdim);
end

function vol = reshape_stack(stack2d, tag)
    npts = size(stack2d, 1);
    nframes = size(stack2d, 2);
    normDim = round(sqrt(npts));
    if normDim^2 ~= npts
        error('get_vCS_TR:BadPlaneSize', '%s: plane size %d is not a perfect square.', tag, npts);
    end
    vol = reshape(stack2d, normDim, normDim, nframes);
end

function [stack2d, nframes] = plane_stack_from_vplanes(Vp, pindex)
    v1 = squeeze(Vp.x(pindex, :, :));
    v2 = squeeze(Vp.y(pindex, :, :));
    v3 = squeeze(Vp.z(pindex, :, :));
    stack2d = single(0.1 * (v1 + v2 + v3));
    if size(stack2d, 1) == 1 && size(stack2d, 2) > 1
        stack2d = stack2d';
    end
    npts = size(stack2d, 1);
    nframes = size(stack2d, 2);
    if nframes < 1
        nframes = 1;
    end
    if npts < 1
        error('get_vCS_TR:EmptyPlane', 'Empty velocity plane at pindex=%s.', mat2str(pindex));
    end
end

function out = resize_plane_stack(vol, imdim)
    [h, w, nf] = size(vol);
    if h == imdim && w == imdim
        out = vol;
        return;
    end
    out = zeros(imdim, imdim, nf, 'like', vol);
    for k = 1:nf
        out(:, :, k) = imresize(vol(:, :, k), [imdim imdim], 'nearest');
    end
end
