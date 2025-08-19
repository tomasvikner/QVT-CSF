% 
function [area_val,diam_val,flowPerHeartCycle_val,maxVel_val,PI_val,RI_val,flowPulsatile_val,...
    velMean_val,VplanesAllx,VplanesAlly,VplanesAllz,r,timeMIPcrossection,segmentFull,...
    vTimeFrameave,MAGcrossection,bnumMeanFlow,bnumStdvFlow,StdvFromMean,Planes, ...
    VplanesCSF, flowCSF, StructCS, CSFSEG, T, CSFROI] ... % new structs for CSF analysis
    = paramMap_params_CSF(filetype,branchList,matrix,timeMIP,vMean,back,...
    BGPCdone,directory,nframes,res,MAG,IDXstart,IDXend,handles, ... 
    StructVols, FileNameFlow) % CSF input 

% PARAMMAP_PARAMS_NEW: Create tangent planes and calculate hemodynamics
%   Based on a sliding threshold segmentation algorithm developed by Carson Hoffman
%   Used by: loadpcvipr.m
%   Dependencies: slidingThreshold.m

%% Tangent Plane Creation
set(handles.TextUpdate,'String','Creating Tangent Planes'); drawnow;
d = 2; % dist. behind/ahead of current pt for tangent plane calc (d=2->5pts)
Tangent_V = zeros(0,3);
for n = 1:max(branchList(:,4))
    branchActual = branchList(branchList(:,4)==n,:);
    dir_temp = zeros(size(branchActual,1),3);
    for i = 1:size(branchActual,1)
        % Extract normal to cross-section
        if i < d+1 %if near 1st endpoint
            dir = (branchActual(i+d,1:3) - branchActual(i,1:3));
        elseif i >= size(branchActual,1)-d %if near 2nd endpoint 
            dir = (branchActual(i,1:3) - branchActual(i-d,1:3));
        else % calculate tangent from d points ahead/behind curr point
            dir = (branchActual(i+d,1:3) - branchActual(i-d,1:3));
        end
        dir_temp(i,:) = dir/norm(dir); %tangent vector with magnitude of 1
    end
    Tangent_V = [Tangent_V;dir_temp]; %#ok<*AGROW> %add all tangents to large list
end

% This will find a normalized vector perpendicular to the tangent vector
[~,idx_max] = max(abs(Tangent_V),[],2); %get max unit along rows
idx_max(idx_max==2) = 1; %flatten to 2D
max_pts = sub2ind(size(Tangent_V),(1:size(Tangent_V,1))',idx_max);
temp = zeros(size(Tangent_V));
temp(max_pts) = 1; %binary matrix of location of max unit vectors
[~,idx_shift] = max(abs(circshift(temp,1,2)),[],2); %rotate (ie x->y,z->x)
shift_pts = sub2ind(size(Tangent_V),(1:size(Tangent_V,1))',idx_shift);
V2 = zeros(size(Tangent_V));
V2(max_pts) = Tangent_V(shift_pts);
V2(shift_pts) = -Tangent_V(max_pts);
N = repmat(sqrt(sum(abs(V2).^2,2)),[1 3]); %repeat vel. magnitude as Nx3
V2 = V2./N;
V3 = cross(Tangent_V,V2); %Third vector that is normalized
% V3,V2,Tangent_V are all orthogonal (i.e. dot( V3(1,:),Tangent_V(1,:) )=0)

%% Interpolate
% Get the full tangent plane for all the points
r = 10; %#ok<NASGU> % size of plane to select from non interpolated data is r*2+1 % PLANE SIZE 
InterpVals = 4; %choose the interpolation between points
Side = r*InterpVals; %creates correct number of points for interpolation
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

dimIM = size(timeMIP);
x = 1:dimIM(1);
y = 1:dimIM(2);
z = 1:dimIM(3);

clear V2 V3 shift_pts temp temp2 idx_max idx_shift Mid x_val y_val z_val 
clear N max_pts d dimIM

% Might use to speed up interpolation time if needed for large matrix sizes
% SE = ones(10,10,10);
% CD_bin_new = imdilate(CD_bin,SE);

%% Interpolation
set(handles.TextUpdate,'String','Interpolating Data');drawnow;
% Get interpolated velocity from 3 directions, multipley w/ tangent vector
v1 = interp3(y,x,z,vMean(:,:,:,1),y_full(:),x_full(:),z_full(:),'linear',0);
v2 = interp3(y,x,z,vMean(:,:,:,2),y_full(:),x_full(:),z_full(:),'linear',0);
v3 = interp3(y,x,z,vMean(:,:,:,3),y_full(:),x_full(:),z_full(:),'linear',0);
v1 = reshape(v1,[length(branchList),(width).^2]);
v2 = reshape(v2,[length(branchList),(width).^2]);
v3 = reshape(v3,[length(branchList),(width).^2]);
temp = zeros([size(v1),3]); % used to hold velocity data information
temp(:,:,1) = bsxfun(@times,v1,Tangent_V(:,1)); % dot product here
temp(:,:,2) = bsxfun(@times,v2,Tangent_V(:,2)); % make veloc. through-plane
temp(:,:,3) = bsxfun(@times,v3,Tangent_V(:,3)); % (mm/s)

% Through-plane SPEED for all points (tangent vector dotted with 3D vel)
vTimeFrameave = sqrt(temp(:,:,1).^2 + temp(:,:,2).^2 + temp(:,:,3).^2); %(mm/s)

% Interpolation for complex difference data
CD_int = interp3(y,x,z,timeMIP,y_full(:),x_full(:),z_full(:),'linear',0);
timeMIPcrossection = reshape(CD_int,[length(branchList),(width).^2]);

% Interpolation for magnitude data
Mag_int = interp3(y,x,z,MAG,y_full(:),x_full(:),z_full(:),'linear',0);
MAGcrossection = reshape(Mag_int,[length(branchList),(width).^2]);

% Add similar data for all fields of volume and CS structs
IntVols = []; 
fns = fields(StructVols);
StructCS = []; 
MaxVals = []; 
for f = 1:numel(fns)
    fn = fns{f};
    vol = StructVols.(fn); 
    intvol = interp3(y,x,z,vol,y_full(:),x_full(:),z_full(:),'linear',0); 
    IntVols.(fn) = intvol; 
    StructCS.(fn) = reshape(intvol,[length(branchList),(width).^2]); % this should be sufficient - csfsegFull not needed
    MaxVals.(fn) = max(intvol(:));
end

clear v1 v2 v3 MAG timeMIP temp CD_int Mag_int vtimeave 

%% In-Plane Segmentation
set(handles.TextUpdate,'String','Performing In-Plane Segmentation'); drawnow;
area_val = zeros(size(Tangent_V,1),1);
diam_val = zeros(size(Tangent_V,1),1);
segmentFull = zeros([length(branchList),(width).^2]);

% Edit for CSF
CSFSEG.bics = segmentFull;
CSFSEG.bcsf = segmentFull;
CSFSEG.cnob = segmentFull; 
CSFSEG.auto = segmentFull;
CSFSEG.madj = segmentFull; 
CSFSEG.full = segmentFull; 

CSFSEG.mthrTrack = zeros(size(Tangent_V,1),1);
CSFSEG.sthrTrack = zeros(size(Tangent_V,1),1);

% Edit threshold list 
[aa, bb] = size(segmentFull);
cc = 3; dd = 3; % temp for evaluating local CSF threshold settings 
CSFROI = [];
CSFROI.bics = zeros(aa, bb, cc, dd);
CSFROI.bcsf = zeros(aa, bb, cc, dd);
CSFROI.cnob = zeros(aa, bb, cc, dd);
CSFROI.flow = zeros(aa, 20, cc, dd);

ecsfc = 0; % count CS where CSF but not blood is empty 
for n = 1:size(Tangent_V,1)
    %%%%%% SLIDING THRESHOLD %%%%%%
    % Get Planes and normalize
    CDSLICE = reshape(timeMIPcrossection(n,:),[(width),(width)]);
    temp = CDSLICE - min(CDSLICE); %shift the minimum to 0
    CDSLICE = temp./max(temp(:)); %now normalize from 0 to 1
    
    velSLICE = reshape(vTimeFrameave(n,:),[(width),(width)]);
    temp = velSLICE - min(velSLICE);
    velSLICE = temp./max(temp(:));
    
    MAGSLICE = reshape(MAGcrossection(n,:),[(width),(width)]);
    temp = MAGSLICE - min(MAGSLICE);
    MAGSLICE = temp./max(temp(:));
    
    weightIMS = [.2 .8 .2]; % Weights = [Mag CD Vel]
    weightIMAGE = (weightIMS(1).*MAGSLICE) + (weightIMS(2).*CDSLICE) + (weightIMS(3).*velSLICE);
    
    step = 0.001;
    UPthresh = 0.8;
    SMf = 90; %smoothing factor
    shiftHM_flag = 0; %do not shift by FWHM
    medFilt_flag = 1; %flag for median filtering of CD image
    [~,segment] = slidingThreshold(weightIMAGE,step,UPthresh,SMf,shiftHM_flag,medFilt_flag);
    areaThresh = round(sum(segment(:)).*0.05); %minimum area to keep
    conn = 6; %connectivity (i.e. 6-pt)
    segment = bwareaopen(segment,areaThresh,conn); %inverse fill holes
    % segment0 = segment; % save initial CSF segment before center segment extraction 

    % Remove all segments not closest to the center
    s = regionprops(segment,'centroid'); %centroids of unique lbls  
    CenterIm = [size(segment,1)/2,size(segment,2)/2]; %loc image center
    Centroids = reshape([s(:).Centroid],[2,length([s(:).Centroid])/2])';
    DisCen = sqrt(sum((Centroids - repmat(CenterIm,[size(Centroids,1),1])).^2,2));
    [~,CenIdx]  = min(DisCen); %find centroid closest to center

    % Fill in the holes and clean up
    [L,Num] = bwlabel(segment); %find centroid index
    LabUse = 1:Num;
    segment = L==LabUse(CenIdx); %cut out other centroids

    % *** LOCAL CSF SEGMENTATION ***
    Slice = [];
    fns = fieldnames(StructCS);
    for fi = 1:numel(fns)
        fn = fns{fi};
        Slice.(fn) = reshape(StructCS.(fn)(n, :), [(width),(width)]); 
        Slice.(fn) = Slice.(fn) / max(Slice.(fn)(:)); % Normalize slice values 
        CSFSEG.(fn)(n, :) = Slice.(fn)(:); % 
    end
 
    STDTHRESH = 0.4; % start with relatively high SD threshold, decrease if segment is not large compared to blood segment 
    BICSTHRESH = 0.2; % start with a low BICS threshold - increase if BICS segment is small compared to blood segment  
    CSFSEG.mthrTrack(n) = BICSTHRESH; 
    [bics, bcsf, cnob, bcube, full, stdThresh, thlist] = segment_pcsf(Slice, BICSTHRESH, STDTHRESH, segment); % send slice struct to local segmentation function 
    CSFSEG.bics(n, :) = bics(:); % blood in CSF scan 
    CSFSEG.bcsf(n, :) = bcsf(:); % binary CSF seg 
    CSFSEG.cnob(n, :) = cnob(:); % csf not blood - from bcsf & dcsf & ~bics 
    CSFSEG.bcube(n, :) = bcube(:); % from T2 CUBE 
    CSFSEG.auto(n, :) = cnob(:); % set auto to cnob 
    CSFSEG.full(n, :) = full(:); % just csf mag > 0 (in FoV) 
    CSFSEG.sthrTrack(n) = stdThresh;

    % temp for evaluating local CSF segmentation
    for b1 = 1:5
        for b3 = 1:5
            CSFROI.bics(n, :, b1, b3) = thlist.bics{b1, b3}(:);
            CSFROI.bcsf(n, :, b1, b3) = thlist.bcsf{b1, b3}(:);
            CSFROI.cnob(n, :, b1, b3) = thlist.cnob{b1, b3}(:);
        end
    end

    if mod(n, 500) == 0 % check that its computing
        disp(['In-plane seg for CS-index: ' num2str(n)]);
    end
    
    % Vessel area measurements
    dArea = (res/10).^2; %pixel size (cm^2)
    area_val(n,1) = sum(segment(:))*dArea*((2*r+1)/(2*r*InterpVals+1))^2;
    segmentFull(n,:) = segment(:); % final blood segment l
    
    % New with ratios of areas. Ratio of smallest inner circle over
    % largest encompassing outer circle (assume circular area). Measure of
    % circularity of vessel (ratio =1 is circle,ratio<1 is irregular shape)
    D = bwdist(~segment); %euclidean distance transform
    Rin = max(D(:)); %distance from center to closest non-zero entry
    [xLoc,yLoc] = find(bwperim(segment)); %get perimeter
    D = pdist2([xLoc,yLoc],[xLoc,yLoc]); %distance b/w perimeter points
    Rout = max(D(:))/2; %radius of largest outer circle
    diam_val(n,1) = Rin^2/Rout^2; %ratio of areas
    diam_val(diam_val==inf) = 0;

end 
clear CDSLICE MAGSLICE temp segment weightIMAGE L LabUse CenIdx Num
disp('In-plane seg. done (paramMap_params_CSF)') % 
disp(['Total CBF & CSF CS count: ' num2str(size(Tangent_V, 1))]);
disp(['Empty-CSF (but not CBF) count: ' num2str(ecsfc)]);

% make sure all csfsegFull segmentations have mcsf > 0
fns = fieldnames(StructCS);
for fi = 1:numel(fns)
    fn = fns{fi};
    CSFSEG.(fn) = CSFSEG.(fn) .* (CSFSEG.mcsf > 0);
end

fieldsToAdd = setdiff(fieldnames(CSFSEG), fieldnames(StructCS));
for i = 1:numel(fieldsToAdd)
    StructCS.(fieldsToAdd{i}) = CSFSEG.(fieldsToAdd{i});
end

%% Extract Time-Resolved Velocities

% Initialize time-resolved hemodynamic parameters 
% Sliding Threshold
flowPulsatile_val = zeros(size(area_val,1),nframes);
maxVelFrame = zeros(size(area_val,1),nframes);
velPulsatile_val = zeros(size(area_val,1),nframes);
bnumMeanFlow = zeros(max(branchList(:,4)),1);
bnumStdvFlow = zeros(max(branchList(:,4)),1);
     
% Initialize time-resolved velocity matrix (not interpolated yet)
VplanesAllx = zeros([length(branchList),(r.*2+1).^2 nframes],'single');
VplanesAlly = zeros([length(branchList),(r.*2+1).^2 nframes],'single');
VplanesAllz = zeros([length(branchList),(r.*2+1).^2 nframes],'single');
VplanesCSF.x = zeros([length(branchList),(r.*2+1).^2 nframes],'single');
VplanesCSF.y = zeros([length(branchList),(r.*2+1).^2 nframes],'single');
VplanesCSF.z = zeros([length(branchList),(r.*2+1).^2 nframes],'single');

% Extract single interp location Idx
ROW = repmat((1:InterpVals:width)',[1 r*2+1]); %replicate up-down
COL = repmat(1:InterpVals*(width):(width)^2,[r*2+1 1])-1; %rep. lf-rt
idCOL = reshape(ROW+COL,[1 numel(ROW)]); %interp query points

% simplify this, new load velocities
vxf = h5read(fullfile(directory,FileNameFlow),'/xcbf');
vyf = h5read(fullfile(directory,FileNameFlow),'/ycbf');
vzf = h5read(fullfile(directory,FileNameFlow),'/zcbf');
startInds = [IDXstart(1),IDXstart(2),IDXstart(3)]; 
endInds = [IDXend(1), IDXend(2), IDXend(3)];
xs = startInds(1); xe = endInds(1); 
ys = startInds(2); ye = endInds(2); 
zs = startInds(3); ze = endInds(3);
vxf = single(vxf(xs:xe, ys:ye, zs:ze, :));
vyf = single(vyf(xs:xe, ys:ye, zs:ze, :));
vzf = single(vzf(xs:xe, ys:ye, zs:ze, :));
% Add the CSF data as well 
cxf = h5read(fullfile(directory,FileNameFlow),'/xcsf'); 
cyf = h5read(fullfile(directory,FileNameFlow),'/ycsf');
czf = h5read(fullfile(directory,FileNameFlow),'/zcsf');
cxf = single(cxf(xs:xe, ys:ye, zs:ze, :));
cyf = single(cyf(xs:xe, ys:ye, zs:ze, :));
czf = single(czf(xs:xe, ys:ye, zs:ze, :));

disp('max of xcbf and xcsf: ')
[max(vxf(:)), max(cxf(:))] %#ok<*NOPRT> 

flowCSF = []; % collect all CSF flow waveforms 
nvoxels = size(Tangent_V,1);
% note actual plane-loop (above was for STD) 
mnz = CSFSEG.mcsf > 0; % if CSF magnitude is zero, we're out of FoV 

% edit additional PCA conditions 
ipca = mnz & CSFSEG.bics > 0;
ipca = ipca & CSFSEG.bcsf > 0;

[aa,bb]=size(mnz);
TR4PCA = zeros(aa, bb, nframes); % Collect for CSF PCA waveform of full FoV
flowCSF.PC1.mean = zeros(nvoxels, nframes);
flowCSF.PC1.median = zeros(nvoxels, nframes);
vp1 = zeros(nvoxels, nframes);
vp2 = zeros(nvoxels, nframes);
vp3 = zeros(nvoxels, nframes);
for j = 1:nframes

    % read blood flow 
    vx = vxf(:, :, :, j);
    vy = vyf(:, :, :, j);
    vz = vzf(:, :, :, j);
    
    % Interpolation of time-resolved velocities
    v1 = interp3(y,x,z,vx,y_full(:),x_full(:),z_full(:),'linear',0);
    v2 = interp3(y,x,z,vy,y_full(:),x_full(:),z_full(:),'linear',0);
    v3 = interp3(y,x,z,vz,y_full(:),x_full(:),z_full(:),'linear',0);
    v1 = reshape(v1,[length(branchList),(width).^2]);
    v2 = reshape(v2,[length(branchList),(width).^2]);
    v3 = reshape(v3,[length(branchList),(width).^2]);

    % check direction deviation between streamline and centerline 
    tp1 = zeros(nvoxels, 1);
    tp2 = zeros(nvoxels, 1);
    tp3 = zeros(nvoxels, 1);
    for k = 1:nvoxels
        vinds = find(segmentFull(k, :)');
        tp1(k) = mean(v1(k, vinds));
        tp2(k) = mean(v2(k, vinds));
        tp3(k) = mean(v3(k, vinds));
    end
    vp1(:, j) = tp1;
    vp2(:, j) = tp2;
    vp3(:, j) = tp3;

    v1 = bsxfun(@times,v1,Tangent_V(:,1)); %dot product here
    v2 = bsxfun(@times,v2,Tangent_V(:,2)); %make velocity through-plane
    v3 = bsxfun(@times,v3,Tangent_V(:,3)); %mm/s)
    VplanesAllx(:,:,j) = v1(:,idCOL); % uninterpolated TR vel. (mm/s)
    VplanesAlly(:,:,j) = v2(:,idCOL);
    VplanesAllz(:,:,j) = v3(:,idCOL);

    % add CSF 
    cx = cxf(:, :, :, j);
    cy = cyf(:, :, :, j);
    cz = czf(:, :, :, j);
    c1 = interp3(y,x,z,cx,y_full(:),x_full(:),z_full(:),'linear',0);
    c2 = interp3(y,x,z,cy,y_full(:),x_full(:),z_full(:),'linear',0);
    c3 = interp3(y,x,z,cz,y_full(:),x_full(:),z_full(:),'linear',0);
    c1 = reshape(c1,[length(branchList),(width).^2]);
    c2 = reshape(c2,[length(branchList),(width).^2]);
    c3 = reshape(c3,[length(branchList),(width).^2]);
    c1 = bsxfun(@times,c1,Tangent_V(:,1)); % dot product here
    c2 = bsxfun(@times,c2,Tangent_V(:,2)); % make velocity through-plane
    c3 = bsxfun(@times,c3,Tangent_V(:,3)); % mm/s)
    VplanesCSF.x(:,:,j) = c1(:,idCOL); % uninterpolated TR vel. (mm/s)
    VplanesCSF.y(:,:,j) = c2(:,idCOL);
    VplanesCSF.z(:,:,j) = c3(:,idCOL);
    if mod(j, 5) == 0
        disp(['Planes set (frame: ' num2str(j) ')']); 
    end

    % Sliding Threshold
    vTimeFrame = segmentFull.*(0.1*(v1 + v2 + v3)); % masked velocity (cm/s)
    vTimeFramerowMean = sum(vTimeFrame,2) ./ sum(vTimeFrame~=0,2); %mean vel
    flowPulsatile_val(:,j) = vTimeFramerowMean.*area_val; %TR flow (ml/s)
    maxVelFrame(:,j) = max(vTimeFrame,[],2); % max vel. each frame (cm/s)
    velPulsatile_val(:,j) = vTimeFramerowMean; % mean vel. each frame (cm/s) 

    % add CSF waveforms from local segments 
    fns = fieldnames(CSFSEG);
    for fi = 1:numel(fns)
        fn = fns{fi};
        csfseg = CSFSEG.(fn) .* mnz; % Exclude out of bounds voxels 
        trCSF.(fn) = csfseg.*(0.1*(c1 + c2 + c3)); % time-resolved CSF using segmentation volume (fn) 
        voxres = 0.5/4; 
        nvox = sum(csfseg, 2);
        acsf = nvox * voxres^2; % mm2
        flowCSF.(fn).median(:, j) = median(trCSF.(fn), 2, 'omitnan') .* acsf;
        flowCSF.(fn).mean(:, j) = mean(trCSF.(fn), 2, 'omitnan') .* acsf;
    end

    % edit for CSF threshold evaluation 
    for b1 = 1:5
        for b3 = 1:5
            csfseg = CSFROI.cnob(:, :, b1, b3); 
            nvox = sum(csfseg, 2); acsf = nvox * voxres^2; %#ok<NASGU> % mm2
            trCSF.temp = csfseg.*(0.1*(c1 + c2 + c3));
            CSFROI.flow(:, j, b1, b3) = mean(trCSF.temp, 2, 'omitnan') .* acsf; 
        end
    end

    temp = sum(flowCSF.cnob.median(:));
    if temp < 1
        disp('Empty CnoB flow')
    end
    npca = sum(ipca(:)); % PCA planeselocity scaling after num CS elems
    TR4PCA(:, :, j) = (ipca/npca) .* (0.1*(c1 + c2 + c3));

end 

% --- MULTIPLANE PCA --- (not used) 
bnums = branchList(:, 4);
noffset = 2; 
MPPCA = 1; 
if MPPCA 
    for n = 1 + noffset : nvoxels - noffset 
        if n == 1
            nrange = 1:3;
        elseif n == 2
            nrange = 2:4; 
        elseif n == nvoxels
            nrange = nvoxels-2:nvoxels; 
        elseif n == nvoxels-1
            nrange = nvoxels-3:nvoxels;
        else
            nrange = n-noffset:n+noffset; 
        end
        brange = bnums(nrange); 
        bmode = mode(brange); 
        binds = find(brange == bmode);
        irange = nrange(binds)'; %#ok<*FNDSB> 
        ninc = numel(irange);
        XPCA = reshape(TR4PCA(irange, :, :), ninc*width^2, nframes);
        [coeff, score] = pca(XPCA', 'NumComponents', 1); %#ok<ASGLU> 
        pc1_waveform = score(:,1); 
        meanWaveform = mean(XPCA, 1, 'omitnan')'; % mean over voxels 
        if corr(pc1_waveform, meanWaveform) < 0
            pc1_waveform = -pc1_waveform; % Flip PC1 if it has negative correlation with the mean waveform
        end
        flowCSF.PC1.median(n, :) = pc1_waveform; 
        flowCSF.PC1.mean(n, :) = pc1_waveform;
    end

% --- SINGLEPLANE PCA ---
% PCA of full FoV (where CSF mag > 0)
else
    for n = 1:nvoxels %#ok<UNRCH> 
        XPCA = reshape(TR4PCA(n, :, :), width^2, nframes);
        [coeff, score] = pca(XPCA', 'NumComponents', 1); %#ok<ASGLU> 
        pc1_waveform = score(:,1); 
        meanWaveform = mean(XPCA, 1, 'omitnan')'; % mean over voxels 
        if corr(pc1_waveform, meanWaveform) < 0
            pc1_waveform = -pc1_waveform; % Flip PC1 if it has negative correlation with the mean waveform
        end
        flowCSF.PC1.median(n, :) = pc1_waveform; 
        flowCSF.PC1.mean(n, :) = pc1_waveform;
    end
end

% calculate tangents based on centerline and velocity
T = [];
T.CL = Tangent_V; 
T.vp1 = mean(vp1, 2);
T.vp2 = mean(vp2, 2);
T.vp3 = mean(vp3, 2);
vnorm = sqrt(T.vp1.^2 + T.vp2.^2 + T.vp3.^2);
T.Vel = [T.vp1, T.vp2, T.vp3]./vnorm; 
disp('T.CL-T.Vel x/y/z correlations: ')
rtang = zeros(3, 1); 
for i = 1:3
    % figure, scatter(abs(T.Vel(:, i)), abs(T.CL(:, i)))
    rtang(i) = corr(abs(T.Vel(:, i)), abs(T.CL(:, i))); %#ok<NOPRT> 
end
mrtang = mean(rtang) %#ok<NOPRT> 
fnTcorr = fullfile(directory, 'FlowTang.txt');
fid = fopen(fnTcorr, 'w');
if fid == -1
    error('Cannot open file: %s', fnTcorr);
end
fprintf(fid, '%.6f\n', mrtang);
fclose(fid);

clear COL ROW idCOL Tangent_V v1 v2 v3 vx vy vz x_full y_full z_full x y z
disp('Time-resolved velocities extracted') % 

%% Compute Hemodynamic Parameters
maxVel_val = max(maxVelFrame,[],2); %max in-plane veloc. for all frames
flowPerHeartCycle_val = sum(flowPulsatile_val,2)./(nframes); %TA flow (ml/s)
velMean_val = sum(velPulsatile_val,2)./(nframes); %TA in-plane velocities
% Pulsatility Index (PI) = (systolic vel - diastolic vel)/(mean vel)
PI_val = abs(max(flowPulsatile_val,[],2) - min(flowPulsatile_val,[],2))./mean(flowPulsatile_val,2);
% Resistivity Index (RI) = (systolic vel - diastolic vel)/(systolic vel)
RI_val = abs(max(flowPulsatile_val,[],2) - min(flowPulsatile_val,[],2))./max(flowPulsatile_val,[],2);

% Mean and standard deviation of flow along all branches
for i=1:max(branchList(:,4))
    idx1 = branchList(:,4)==i; %find all points along branch
    bnumMeanFlow(i) = mean(flowPerHeartCycle_val(idx1)); %mean TA flow
    bnumStdvFlow(i) = std(flowPerHeartCycle_val(idx1)); %stdv TA flow
end 

% Get coefficient of variation (stdv from mean) for all points along branch
% Looks at local stdv and mean (window width of 5).
StdvFromMean = flowPerHeartCycle_val;
for n = 1:max(branchList(:,4))
    IDbranch = find(branchList(:,4)== n); %extract points for branch n
    % Calculate near branch start
    StdvFromMean(IDbranch(1)) = std(flowPerHeartCycle_val(IDbranch(1:3))) ./ abs(mean(flowPerHeartCycle_val(IDbranch(1:3))));
    StdvFromMean(IDbranch(2)) = std(flowPerHeartCycle_val(IDbranch(1:4))) ./ abs(mean(flowPerHeartCycle_val(IDbranch(1:4))));
    % Calculate for middle of branch (window width of 5)
    for m = 1:numel(IDbranch)-4
        StdvFromMean(IDbranch(m+2)) = std(flowPerHeartCycle_val(IDbranch(m:m+4)))./abs(mean(flowPerHeartCycle_val(IDbranch(m:m+4))));
    end
    % Calculate near branch end
    StdvFromMean(IDbranch(end-1)) = std(flowPerHeartCycle_val(IDbranch(end-3:end)))./abs(mean(flowPerHeartCycle_val(IDbranch(end-3:end))));
    StdvFromMean(IDbranch(end)) = std(flowPerHeartCycle_val(IDbranch(end-2:end)))./abs(mean(flowPerHeartCycle_val(IDbranch(end-2:end))));
end
StdvFromMean = StdvFromMean - min(StdvFromMean(:)); %shift the minimum to 0
StdvFromMean = StdvFromMean./max(StdvFromMean(:)); %normalize range 0-1
end
