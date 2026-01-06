function [VcrossTR, CcrossTR] = get_vCS_TR(Vplanes, pindex)

    % Orig CS data 
    v1 = squeeze(Vplanes.CBF.x(pindex,:,:));
    v2 = squeeze(Vplanes.CBF.y(pindex,:,:));
    v3 = squeeze(Vplanes.CBF.z(pindex,:,:));
    VcrossTR = 0.1*(v1 + v2 + v3); 

    normDim = sqrt(size(VcrossTR,1));
    VcrossTR = reshape(VcrossTR,normDim,normDim,nframes);
    VcrossTR = imresize(VcrossTR,[imdim imdim],'nearest');

    % CSF CS data
    c1 = squeeze(Vplanes.CSF.x(pindex,:,:));
    c2 = squeeze(Vplanes.CSF.y(pindex,:,:));
    c3 = squeeze(Vplanes.CSF.z(pindex,:,:));
    CcrossTR = 0.1*(c1 + c2 + c3); 

    normDim = sqrt(size(CcrossTR,1));
    CcrossTR = reshape(CcrossTR,normDim,normDim,nframes);
    CcrossTR = imresize(CcrossTR,[imdim imdim],'nearest');

end