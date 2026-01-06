function flowCSF = compute_PCA_CSF(flowCSF, TR4PCA, branchList, ...
                                  width, nframes, nvoxels, ...
                                  MULTIPLANE, noffset)
% COMPUTE_PCA_CSF
% Computes PC1 waveforms for CSF flow using either:
%   - Multiplane PCA (branch-consistent neighborhood)
%   - Single-plane PCA (per-voxel plane)
%
% OUTPUT:
%   flowCSF.PC1.mean   [nvoxels x nframes]
%   flowCSF.PC1.median[nvoxels x nframes]

% Preallocate
flowCSF.PC1.mean   = nan(nvoxels, nframes);
flowCSF.PC1.median = nan(nvoxels, nframes);

bnums = branchList(:,4);

for n = 1:nvoxels

    % ================= MULTIPLANE =================
    if MULTIPLANE

        if n <= noffset
            nrange = 1 : min(2*noffset+1, nvoxels);
        elseif n > nvoxels-noffset
            nrange = max(1, nvoxels-2*noffset) : nvoxels;
        else
            nrange = n-noffset : n+noffset;
        end

        brange = bnums(nrange);
        bmode  = mode(brange);
        irange = nrange(brange == bmode)';

        XPCA = reshape(TR4PCA(irange,:,:), ...
                        numel(irange)*width^2, nframes);

    % ================= SINGLEPLANE =================
    else
        XPCA = reshape(TR4PCA(n,:,:), width^2, nframes);
    end

    % ================= PCA =================
    [~, score] = pca(XPCA', 'NumComponents', 1);
    pc1 = score(:,1);

    meanWaveform = mean(XPCA, 1, 'omitnan')';
    if corr(pc1, meanWaveform) < 0
        pc1 = -pc1;
    end

    flowCSF.PC1.mean(n,:)   = pc1;
    flowCSF.PC1.median(n,:) = pc1;
end

end
