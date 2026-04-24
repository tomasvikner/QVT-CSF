function [Vplanes, v1, v2, v3, c1, c2, c3] = calc_vplanes( ...
    vxf, vyf, vzf, ...
    cxf, cyf, czf, ...
    x, y, z, x_full, y_full, z_full, ...
    Tangent_V, branchList, ...
    width, r, InterpVals, nframes)
xg = 1:x; yg = 1:y; zg = 1:z;
if ndims(vxf) >= 4
    nframesData = min([size(vxf,4), size(vyf,4), size(vzf,4), ...
        size(cxf,4), size(cyf,4), size(czf,4)]);
else
    nframesData = 1;
end
nframes = max(1, min(nframes, nframesData));

% ==========================================================
% INTERP COLUMN INDICES
% ==========================================================
ROW   = repmat((1:InterpVals:width)',[1 2*r+1]);
COL   = repmat(1:InterpVals*width:width^2,[2*r+1 1]) - 1;
idCOL = reshape(ROW + COL, 1, []);

nvox  = size(branchList,1);
npts  = numel(idCOL);
npix  = width^2;

% ==========================================================
% OUTPUT
% ==========================================================
Vplanes.CBF.x = zeros(nvox,npts,nframes,'single');
Vplanes.CBF.y = zeros(nvox,npts,nframes,'single');
Vplanes.CBF.z = zeros(nvox,npts,nframes,'single');
Vplanes.CSF.x = zeros(nvox,npts,nframes,'single');
Vplanes.CSF.y = zeros(nvox,npts,nframes,'single');
Vplanes.CSF.z = zeros(nvox,npts,nframes,'single');

v1 = zeros(nvox,npix,nframes,'single');
v2 = v1; v3 = v1;
c1 = v1; c2 = v1; c3 = v1;

% ==========================================================
% FRAME LOOP
% ==========================================================
for j = 1:nframes

    % -------- CBF --------
    v1j = interp3(yg,xg,zg,vxf(:,:,:,j),y_full(:),x_full(:),z_full(:),'linear',0);
    v2j = interp3(yg,xg,zg,vyf(:,:,:,j),y_full(:),x_full(:),z_full(:),'linear',0);
    v3j = interp3(yg,xg,zg,vzf(:,:,:,j),y_full(:),x_full(:),z_full(:),'linear',0);

    v1j = reshape(v1j,nvox,npix) .* Tangent_V(:,1);
    v2j = reshape(v2j,nvox,npix) .* Tangent_V(:,2);
    v3j = reshape(v3j,nvox,npix) .* Tangent_V(:,3);

    v1(:,:,j) = v1j;  v2(:,:,j) = v2j;  v3(:,:,j) = v3j;

    Vplanes.CBF.x(:,:,j) = v1j(:,idCOL);
    Vplanes.CBF.y(:,:,j) = v2j(:,idCOL);
    Vplanes.CBF.z(:,:,j) = v3j(:,idCOL);

    % -------- CSF --------
    c1j = interp3(yg,xg,zg,cxf(:,:,:,j),y_full(:),x_full(:),z_full(:),'linear',0);
    c2j = interp3(yg,xg,zg,cyf(:,:,:,j),y_full(:),x_full(:),z_full(:),'linear',0);
    c3j = interp3(yg,xg,zg,czf(:,:,:,j),y_full(:),x_full(:),z_full(:),'linear',0);

    c1j = reshape(c1j,nvox,npix) .* Tangent_V(:,1);
    c2j = reshape(c2j,nvox,npix) .* Tangent_V(:,2);
    c3j = reshape(c3j,nvox,npix) .* Tangent_V(:,3);

    c1(:,:,j) = c1j;  c2(:,:,j) = c2j;  c3(:,:,j) = c3j;

    Vplanes.CSF.x(:,:,j) = c1j(:,idCOL);
    Vplanes.CSF.y(:,:,j) = c2j(:,idCOL);
    Vplanes.CSF.z(:,:,j) = c3j(:,idCOL);

    if mod(j,5)==0
        disp(['Planes set (frame: ' num2str(j) ')']);
    end

end
end
