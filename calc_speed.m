function vTimeFrameave = calc_speed(x, y, z, vMean, x_full, y_full, z_full, Tangent_V)

% Get interpolated velocity from 3 directions, multipley w/ tangent vector
xg = 1:x; yg = 1:y; zg = 1:z;
v1 = interp3(yg,xg,zg,vMean(:,:,:,1),y_full(:),x_full(:),z_full(:),'linear',0);
v2 = interp3(yg,xg,zg,vMean(:,:,:,2),y_full(:),x_full(:),z_full(:),'linear',0);
v3 = interp3(yg,xg,zg,vMean(:,:,:,3),y_full(:),x_full(:),z_full(:),'linear',0);
npts = size(x_full,1);
npx = size(x_full,2);
v1 = reshape(v1,[npts, npx]);
v2 = reshape(v2,[npts, npx]);
v3 = reshape(v3,[npts, npx]);
temp = zeros([size(v1),3]); % used to hold velocity data information
temp(:,:,1) = bsxfun(@times,v1,Tangent_V(:,1)); % dot product here
temp(:,:,2) = bsxfun(@times,v2,Tangent_V(:,2)); % make veloc. through-plane
temp(:,:,3) = bsxfun(@times,v3,Tangent_V(:,3)); % (mm/s)

% Through-plane SPEED for all points (tangent vector dotted with 3D vel)
vTimeFrameave = sqrt(temp(:,:,1).^2 + temp(:,:,2).^2 + temp(:,:,3).^2); %(mm/s)

end