% first row (vals) to struct CLVALS
function [CLVALS, Planes, Vplanes, CBFSEG, CSFSEG, StructCS, CSFROI, MD, flowCSF] ...
    = paramMap_params_CSF(MD, branchList, CD, vMean, MAG, IDX, handles, VOLS, vel4dResampled)

% PARAMMAP_PARAMS_CSF: Modified 2025 for CBF-CSF coupling (Tomas Vikner)
% Structs to pass between functions:
% CLVALS: hemodynamic/CSF parameters along centerline
% Planes:
% Vplanes: CBF and CSF velocities through all planes
% CSFSEG: Local CSF seg across CL (npoints, width^2); multiple options
% CBFSEG: Local CBF seg across CL (prev. segmentFull); CD-based only
% StructCS: TEMP: Similar to CBFSEG?
% CSFROI: CSF flow data as a function of CNOB seg settings (b1, b3)
% MD: Metadata, including flow file name, directory name, and other info
% IDX: IDX.start and IDXend. for bounding box
% VOLS: 3D volumes in addition to CD/MAG from loadHDF5.m
% vel4dResampled: optional struct with CBF+CSF resampled 4-D fields vxf,vyf,vzf,cxf,cyf,czf
%   (from loadHDF5 only; not part of MD). Omit or pass [] to load/crop from HDF5 instead.

if nargin < 9 || isempty(vel4dResampled)
    useResampled4d = false;
else
    useResampled4d = isstruct(vel4dResampled) && isfield(vel4dResampled, 'vxf');
end

% Tangent Plane Creation
if ~isfield(MD, 'RECONRES') || isempty(MD.RECONRES)
    MD.RECONRES = 0.75; % mm isotropic standard reconstruction resolution
end
if ~isfield(MD, 'iframes') || isempty(MD.iframes)
    MD.iframes = 1000; % standard high-res phase grid (interpCoupling; set_metadata also sets this)
end

set(handles.TextUpdate,'String','Creating Tangent Planes'); drawnow;
[Tangent_V, V2, V3] = calc_tangent(branchList); %#ok<*ASGLU>

% Interpolate and get the full tangent plane for all the points
PLANESIZE = 10; % size of plane to select from non interpolated data is r*2+1 % PLANE SIZE 
RECONRES = MD.RECONRES;
InterpVals = 4; % choose the interpolation between points
set(handles.TextUpdate,'String','Interpolating Data'); drawnow;
[x_full, y_full, z_full, Planes] = calc_planes(PLANESIZE, InterpVals, branchList, V2, V3); 
[x, y, z] = size(CD);
vTimeFrameave = calc_speed(x, y, z, vMean, x_full, y_full, z_full, Tangent_V); 
xg = 1:x; yg = 1:y; zg = 1:z;
width = round(sqrt(size(x_full,2))); % plane width in pixels
r = PLANESIZE; % legacy variable name used downstream

% Interpolation for complex difference and magnitude data
StructCS = []; 
CD_int = interp3(yg,xg,zg,CD,y_full(:),x_full(:),z_full(:),'linear',0);
Mag_int = interp3(yg,xg,zg,MAG,y_full(:),x_full(:),z_full(:),'linear',0);
StructCS.MAG = reshape(Mag_int,[length(branchList),(width).^2]);
StructCS.CD = reshape(CD_int,[length(branchList),(width).^2]);

% Add similar data for all fields of volume and CS structs
IntVols = []; 
fns = fields(VOLS);
MaxVals = []; 
for f = 1:numel(fns)
    fn = fns{f};
    vol = VOLS.(fn); 
    intvol = interp3(yg,xg,zg,vol,y_full(:),x_full(:),z_full(:),'linear',0); 
    IntVols.(fn) = intvol; 
    StructCS.(fn) = reshape(intvol,[length(branchList),(width).^2]); % this should be sufficient - csfsegFull not needed
    MaxVals.(fn) = max(intvol(:));
end
if ~isfield(StructCS, 'cube')
    StructCS.cube = ones(size(StructCS.mcsf), 'like', StructCS.mcsf);
end

% In-Plane Segmentation
CLVALS = []; 
set(handles.TextUpdate,'String','Performing In-Plane Segmentation'); drawnow;
area_val = zeros(size(Tangent_V,1),1);
diam_val = zeros(size(Tangent_V,1),1);

% Allocate for local CBF and CSF seg
base = zeros(length(branchList), width^2);
CBFSEG = base;
CSFSEG = struct( ...
    'bics', base, ...
    'bcsf', base, ...
    'cnob', base, ...
    'auto', base, ...
    'madj', base, ...
    'full', base );

% Magnitude and vel SD thresholds track 
CSFSEG.mthrTrack = zeros(size(Tangent_V,1),1);
CSFSEG.sthrTrack = zeros(size(Tangent_V,1),1);

% Edit threshold list 
[aa, bb] = size(CBFSEG);
cc = 3; dd = 3; 
CSFROI = []; % TEMP: evaluate local CSF threshold settings (2025 paper)
CSFROI.bics = zeros(aa, bb, cc, dd);
CSFROI.bcsf = zeros(aa, bb, cc, dd);
CSFROI.cnob = zeros(aa, bb, cc, dd);

ecsfc = 0; % temp: debug/count CS where CSF (but not blood) is empty 
for n = 1:size(Tangent_V,1)

    % *** LOCAL CBF SEGMENTATION (slidingThreshold) ***
    norm01 = @(x) (x - min(x(:))) ./ max(eps, max(x(:)) - min(x(:)));
    CDSLICE  = norm01(reshape(StructCS.CD(n,:),  width, width));
    velSLICE = norm01(reshape(vTimeFrameave(n,:), width, width));
    MAGSLICE = norm01(reshape(StructCS.MAG(n,:), width, width));
    
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
    segment = bwareaopen(segment,areaThresh,conn); % inverse fill holes

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
    
    % Vessel area measurements
    dArea = (MD.RECONRES/10).^2; % pixel size (cm^2), based on recon resolution
    area_val(n,1) = sum(segment(:))*dArea*((2*r+1)/(2*r*InterpVals+1))^2;
    CBFSEG(n,:) = segment(:); % final blood segment l
    diam_val(n,1) = calc_diam(segment); % TEMP note: is this really diam? 
    % *** END OF LOCAL CBF SEGMENTATION ***

    % *** LOCAL CSF SEGMENTATION (QVT-CSF only) ***
    Slice = [];
    fns = fieldnames(StructCS);
    for fi = 1:numel(fns)
        fn = fns{fi};
        Slice.(fn) = reshape(StructCS.(fn)(n, :), [(width),(width)]); 
        Slice.(fn) = Slice.(fn) / max(Slice.(fn)(:)); % Normalize slice values 
        CSFSEG.(fn)(n, :) = Slice.(fn)(:); % 
    end
    if ~isfield(Slice, 'cube')
        Slice.cube = ones(width, width, 'like', Slice.mcsf);
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

    % Temp for evaluating local CSF segmentation
    for b1 = 1:5
        for b3 = 1:5
            CSFROI.bics(n, :, b1, b3) = thlist.bics{b1, b3}(:);
            CSFROI.bcsf(n, :, b1, b3) = thlist.bcsf{b1, b3}(:);
            CSFROI.cnob(n, :, b1, b3) = thlist.cnob{b1, b3}(:);
        end
    end
    % *** END OF LOCAL CSF SEGMENTATION ***

    % Check for running 
    if mod(n,500)==0, disp(['In-plane seg for CS ind: ' num2str(n)]); end

end 
disp('In-plane seg. done (paramMap_params_CSF)') 
disp(['Total CBF and CSF CS-count: ' num2str(size(Tangent_V, 1))]);
if ecsfc > 0, disp(['Empty CSF (but not CBF) count: ' num2str(ecsfc)]); end

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

% ==========================================================
% LOAD + CROP
% ==========================================================
xs = IDX.start(1); xe = IDX.end(1);
ys = IDX.start(2); ye = IDX.end(2);
zs = IDX.start(3); ze = IDX.end(3);

if useResampled4d
    % Resampled 4-D from loadHDF5 (isotropic grid + temporal resample); not read from MD
    vxf = vel4dResampled.vxf;
    vyf = vel4dResampled.vyf;
    vzf = vel4dResampled.vzf;
    cxf = vel4dResampled.cxf;
    cyf = vel4dResampled.cyf;
    czf = vel4dResampled.czf;
    if ndims(vxf) >= 4
        nframesData = min([size(vxf, 4), size(vyf, 4), size(vzf, 4), ...
            size(cxf, 4), size(cyf, 4), size(czf, 4)]);
    else
        nframesData = 1;
    end
else
    read4d = @(name) single(h5read(fullfile(MD.directory, MD.infile), name));

    % CBF (blood): /bphx, /bphy, /bphz
    vxf = read4d('/bphx'); vxf = vxf(xs:xe, ys:ye, zs:ze, :);
    vyf = read4d('/bphy'); vyf = vyf(xs:xe, ys:ye, zs:ze, :);
    vzf = read4d('/bphz'); vzf = vzf(xs:xe, ys:ye, zs:ze, :); % CBF

    % CSF: /cphx, /cphy, /cphz
    cxf = read4d('/cphx'); cxf = cxf(xs:xe, ys:ye, zs:ze, :);
    cyf = read4d('/cphy'); cyf = cyf(xs:xe, ys:ye, zs:ze, :);
    czf = read4d('/cphz'); czf = czf(xs:xe, ys:ye, zs:ze, :); % CSF
    if ndims(vxf) >= 4
        nframesData = min([size(vxf, 4), size(vyf, 4), size(vzf, 4), ...
            size(cxf, 4), size(cyf, 4), size(czf, 4)]);
    else
        nframesData = 1;
    end
end
if ~isfield(MD, 'nframes') || isempty(MD.nframes)
    MD.nframes = nframesData;
else
    MD.nframes = min(MD.nframes, nframesData);
end
if MD.nframes < 1
    MD.nframes = 1;
end

CSFROI.flow = zeros(size(CSFROI.cnob, 1), MD.nframes, cc, dd);

mnz = StructCS.mcsf > 0; % in-FOV CSF magnitude mask
[aa,bb]=size(mnz); 
ipca = (mnz & CSFSEG.bics > 0) & (CSFSEG.bcsf > 0); 
TR4PCA = zeros(aa, bb, MD.nframes); % Collect for CSF-PCA waveform

% Calculate CBF and CSF velocity planes to struct Vplanes
[Vplanes, v1, v2, v3, c1, c2, c3] = calc_vplanes(vxf, vyf, vzf, ...
    cxf, cyf, czf, x, y, z, x_full, y_full, z_full, Tangent_V, ...
    branchList, width, PLANESIZE, InterpVals, MD.nframes);

% ==========================================================
% Allocate CLVALS + helpers
% ==========================================================
nvox   = size(area_val,1);
nfrm   = MD.nframes;
nb     = max(branchList(:,4));

CLVALS.flowTR   = zeros(nvox,nfrm);
CLVALS.maxvel  = zeros(nvox,1);
CLVALS.meanvel = zeros(nvox,1);
CLVALS.flowTA  = zeros(nvox,1);
CLVALS.PI      = zeros(nvox,1);
CLVALS.RI      = zeros(nvox,1);

maxVelFrame = zeros(nvox,nfrm);
bnumMeanFlow = zeros(nb,1);
bnumStdvFlow = zeros(nb,1);

voxres = RECONRES/InterpVals;
fns = fieldnames(CSFSEG);
acsfByFn = struct();
for fi = 1:numel(fns)
    fn = fns{fi};
    acsfByFn.(fn) = sum(CSFSEG.(fn) .* mnz, 2) * voxres^2;
end
acsfROI = zeros(size(CSFROI.cnob,1), 5, 5);
for b1 = 1:5
    for b3 = 1:5
        acsfROI(:,b1,b3) = sum(CSFROI.cnob(:,:,b1,b3),2) * voxres^2;
    end
end

% ==========================================================
% Time-resolved loop
% ==========================================================
for j = 1:nfrm

    % ---------------- CBF FLOW ----------------
    vTime = CBFSEG .* (0.1*(v1(:,:,j) + v2(:,:,j) + v3(:,:,j)));
    vmean = sum(vTime,2) ./ sum(vTime~=0,2);
    vmean = squeeze(vmean);

    CLVALS.flowTR(:,j) = vmean .* area_val;
    maxVelFrame(:,j)   = squeeze(max(vTime,[],2));

    % ---------------- CSF SEGMENTS ----------------
    for fi = 1:numel(fns)
        fn = fns{fi};
        csfseg = CSFSEG.(fn) .* mnz;
        trCSF  = csfseg .* (0.1*(c1(:,:,j) + c2(:,:,j) + c3(:,:,j)));
        acsf   = acsfByFn.(fn);

        CLVALS.flowCSF.(fn).median(:,j) = squeeze(median(trCSF,2,'omitnan')) .* acsf;
        CLVALS.flowCSF.(fn).mean(:,j)   = squeeze(mean(trCSF,2,'omitnan'))   .* acsf;
    end

    % ---------------- ROI THRESHOLD ----------------
    for b1 = 1:5
        for b3 = 1:5
            csfseg = CSFROI.cnob(:,:,b1,b3);
            acsf   = acsfROI(:,b1,b3);
            trCSF  = csfseg .* (0.1*(c1(:,:,j) + c2(:,:,j) + c3(:,:,j)));
            CSFROI.flow(:,j,b1,b3) = squeeze(mean(trCSF,2,'omitnan')) .* acsf;
        end
    end

    % ---------------- PCA INPUT ----------------
    TR4PCA(:,:,j) = (ipca/sum(ipca(:))) .* (0.1*(c1(:,:,j) + c2(:,:,j) + c3(:,:,j)));

    if mod(j,5)==0
        disp(['Time-resolved velocities extracted (frame ' num2str(j) ')']);
    end
end

% ==========================================================
% Derived metrics (single-pass)
% ==========================================================
CLVALS.maxvel  = max(maxVelFrame,[],2);
CLVALS.flowTA  = mean(CLVALS.flowTR,2);
CLVALS.meanvel = mean(CLVALS.flowTR ./ area_val, 2);

CLVALS.PI = abs(max(CLVALS.flowTR,[],2) - min(CLVALS.flowTR,[],2)) ...
            ./ mean(CLVALS.flowTR,2);

CLVALS.RI = abs(max(CLVALS.flowTR,[],2) - min(CLVALS.flowTR,[],2)) ...
            ./ max(CLVALS.flowTR,[],2);

% ==========================================================
% Branch-wise statistics
% ==========================================================
for b = 1:nb
    idx = branchList(:,4)==b;
    bnumMeanFlow(b) = mean(CLVALS.flowTA(idx));
    bnumStdvFlow(b) = std(CLVALS.flowTA(idx));
end

% ==========================================================
% Local flow consistency (StdvFromMean)
% ==========================================================
StdvFromMean = CLVALS.flowTA;
for b = 1:nb
    id = find(branchList(:,4)==b);
    for k = 1:numel(id)
        w = max(1,k-2):min(numel(id),k+2);
        StdvFromMean(id(k)) = std(CLVALS.flowTA(id(w))) ...
                              ./ abs(mean(CLVALS.flowTA(id(w))));
    end
end

StdvFromMean = (StdvFromMean - min(StdvFromMean)) ./ max(StdvFromMean);

% ==========================================================
% Store remaining CL / MD values
% ==========================================================
CLVALS.area     = area_val;
CLVALS.diam     = diam_val;
CLVALS.StdvFromMean = StdvFromMean;

MD.bnumMeanFlow = bnumMeanFlow;
MD.bnumStdvFlow = bnumStdvFlow;
MD.StdvFromMean = StdvFromMean;
MD.PLANESIZE    = PLANESIZE;
MD.RECONRES     = RECONRES; 

% Waveform struct for datacursor / calc_waveforms (PC1 + per-segment mean/median)
flowCSF = CLVALS.flowCSF;
nvoxPCA = size(branchList,1);
try
    flowCSF = compute_PCA_CSF(flowCSF, TR4PCA, branchList, width, MD.nframes, nvoxPCA, false, 0);
catch ME
    warning('paramMap_params_CSF:PCAFailed', 'compute_PCA_CSF failed: %s', ME.message);
    flowCSF.PC1.mean = nan(nvoxPCA, MD.nframes);
    flowCSF.PC1.median = flowCSF.PC1.mean;
end

end
