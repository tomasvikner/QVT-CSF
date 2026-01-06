function [x_full, y_full, z_full, Planes] = calc_planes(PLANESIZE, InterpVals, branchList, V2)

Side = PLANESIZE*InterpVals; %creates correct number of points for interpolation
width = Side.*2+1; %width of plane in pixels
Mid = zeros(length(branchList),1);

% Find x values on line
temp = repmat(V2(:,1)./InterpVals,[1 Side]);
temp = cumsum(temp,2); %runs from 0 to +(r*interpVals) by unit dist/interp
temp2 = -fliplr(temp); %runs from -(r*interpVals) to 0 by unit dist/interp
x_val = [temp2 Mid temp]; %combine temps--size = N x (r*interpVals*2)+1
x_val = bsxfun(@plus,x_val,branchList(:,1)); %pointwise addition
x_val = reshape(x_val,[numel(x_val) 1]); %stretch into vector

% Find y values on line
temp = repmat(V2(:,2)./InterpVals,[1 Side]);
temp = cumsum(temp,2);
temp2 = -fliplr(temp);
y_val = [temp2 Mid temp];
y_val = bsxfun(@plus,y_val,branchList(:,2));
y_val = reshape(y_val,[numel(y_val) 1]);

% Find z values on the line
temp = repmat(V2(:,3)./InterpVals,[1 Side]);
temp = cumsum(temp,2);
temp2 = -fliplr(temp);
z_val = [temp2 Mid temp];
z_val = bsxfun(@plus,z_val,branchList(:,3));
z_val = reshape(z_val,[numel(z_val) 1]);

% At this point x,y,z values have created a tangent line perpendicular to
% the normal vector for all centerline points.
% Now, we begin filling out the other perpendicular line to create a plane.

% Find x values on plane
Mid = zeros(length(branchList)*(width),1);
temp = repmat(V3(:,1)./InterpVals,[width Side]);
temp = cumsum(temp,2);
temp2 = -fliplr(temp);
x_full = [temp2 Mid temp];
x_full = bsxfun(@plus,x_full,x_val);
x_full = reshape(x_full,[length(branchList)*(width).^2,1]);

% Find y values on plane
temp = repmat(V3(:,2)./InterpVals,[(width) Side]);
temp = cumsum(temp,2);
temp2 = -fliplr(temp);
y_full = [temp2 Mid temp];
y_full = bsxfun(@plus,y_full,y_val);
y_full = reshape(y_full,[length(branchList)*(width).^2,1]);

% Find z values on plane
temp = repmat(V3(:,3)./InterpVals,[(width) Side]);
temp = cumsum(temp,2);
temp2 = -fliplr(temp);
z_full = [temp2 Mid temp];
z_full = bsxfun(@plus,z_full,z_val);
z_full = reshape(z_full,[length(branchList)*(width).^2,1]);

% Typecast to single and reshape
x_full = reshape(single(x_full),[length(branchList),(width).^2]);
y_full = reshape(single(y_full),[length(branchList),(width).^2]);
z_full = reshape(single(z_full),[length(branchList),(width).^2]);

% Get corners of UNINTERPOLATED planes
Planes = zeros(size(branchList,1),4,3);
Planes(:,:,1) = [x_full(:,1),x_full(:,width-InterpVals),x_full(:,end),x_full(:,end-width+1)];
Planes(:,:,2) = [y_full(:,1),y_full(:,width-InterpVals),y_full(:,end),y_full(:,end-width+1)];
Planes(:,:,3) = [z_full(:,1),z_full(:,width-InterpVals),z_full(:,end),z_full(:,end-width+1)];

end