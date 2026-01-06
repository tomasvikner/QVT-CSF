% function [nframes,matrix,res,timeres,VENC,area_val,diam_val,flowPerHeartCycle_val, ...
%     maxVel_val,PI_val,RI_val,flowPulsatile_val,velMean_val, ...
%     VplanesAllx,VplanesAlly,VplanesAllz,Planes,branchList,segment,r, ...
%     timeMIPcrossection,segmentFull,vTimeFrameave,MAGcrossection, imageData, ...
%     bnumMeanFlow,bnumStdvFlow,StdvFromMean, ... 
%     VplanesCSF, flowCSF, StructCS, CSFSEG, T, CSFROI] = loadHDF5(METADATA.directory,handles,METADATA.infile)

% TODO: fix return
function [CLVALS, Planes, Vplanes, CBFSEG, CSFSEG, StructCS, CSFROI, METADATA, imageData] = loadHDF5(METADATA, handles)

% LOADHDF5: loadhdf5 reads in PCVIPR data saved in h5
%   Used by: paramMap.m
%   Dependencies: background_phase_correction.m, evaluate_poly.m, calc_angio.m,
%   feature_extraction.m, paramMap_params_new.m, slidingThreshold.m

%% Read HDF5 (update for CBF/CSF combined .h5-files) 
set(handles.TextUpdate,'String','Loading .HDF5 Data'); 
drawnow;
 
% Read PCVIPR header 
headerFile = [METADATA.directory '/pcvipr_header.txt']; 
PCVIPR_HEADER = read_pcvipr_header(headerFile); 

% add all anatomical information to a volume struct 
StructVols = []; 
StructVols.mcbf = h5read(fullfile(METADATA.directory, METADATA.infile),'/mcbf');
StructVols.mcsf = h5read(fullfile(METADATA.directory, METADATA.infile),'/mcsf');
StructVols.scsf = h5read(fullfile(METADATA.directory, METADATA.infile),'/scsf');
StructVols.scbf = h5read(fullfile(METADATA.directory, METADATA.infile),'/scbf');
StructVols.cube = h5read(fullfile(METADATA.directory, METADATA.infile),'/cube'); 
StructVols.rage = h5read(fullfile(METADATA.directory, METADATA.infile),'/rage');
StructVols.pvas = h5read(fullfile(METADATA.directory, METADATA.infile),'/pvasc');
StructVols.gcsf = h5read(fullfile(METADATA.directory, METADATA.infile),'/gcsf');

% Calculate start and end inds for cropping 3D and 4D data
MAG = single(StructVols.mcbf);
IDX = calc_crop(MAG);

% Crop MAG and average velocity 
MAG = MAG(IDX.start(1):IDX.end(1),IDX.start(2):IDX.end(2),IDX.start(3):IDX.end(3));
newDIM = size(MAG);
vMean = load_vMean(newDIM, METADATA); % Load cropped mean vel. based on dims and metadata (METADATA) 

% Crop all other 3D volumes
fns = fields(StructVols);
for j = 1:numel(fns) 
    fn = fns{j};
    StructVols.(fn) = StructVols.(fn)(IDX.start(1):IDX.end(1),IDX.start(2):IDX.end(2),IDX.start(3):IDX.end(3));
end

% Manual Background Phase Correction (if necessary)
% NOTE: Currently, this is not applied to TR-CSF, needs fix if used
BGPCdone = 1; % if 1, assume BG corr was done in recon; 0 otherwise
back = zeros(size(vMean),'single');
if ~BGPCdone
    set(handles.TextUpdate,'String','Phase Correction with Polynomial'); drawnow; %#ok<*UNRCH>
    
    [poly_fitx,poly_fity,poly_fitz] = background_phase_correction(MAG,vMean(:,:,:,1),vMean(:,:,:,2),vMean(:,:,:,3));
    disp('Correcting data with polynomial');
    xrange = single(linspace(-1,1,size(MAG,1)));
    yrange = single(linspace(-1,1,size(MAG,2)));
    zrange = single(linspace(-1,1,size(MAG,3)));
    [Y,X,Z] = meshgrid(yrange,xrange,zrange);
    
    % Get poly data and correct average velocity for x,y,z dimensions
    back(:,:,:,1) = single(evaluate_poly(X,Y,Z,poly_fitx));
    back(:,:,:,2) = single(evaluate_poly(X,Y,Z,poly_fity));
    back(:,:,:,3) = single(evaluate_poly(X,Y,Z,poly_fitz));
    vMean = vMean - back;
end

% Calculate complex difference angiogram for visualization.
set(handles.TextUpdate,'String','Creating Angiogram'); drawnow;
timeMIP = calc_angio(MAG, vMean, PCVIPR_HEADER.VENC); % approximate CD.dat

% Find optimum global threshold 
step = 0.001; % step size for sliding threshold
UPthresh = 0.8; % max upper threshold when creating Sval curvature plot
SMf = 10; % smoothing factor
shiftHM_flag = 1; % flag to shift max curvature by FWHM
medFilt_flag = 1; % flag for median filtering of CD image
[~,segment] = slidingThreshold(timeMIP,step,UPthresh,SMf,shiftHM_flag,medFilt_flag);
areaThresh = round(sum(segment(:)).*0.01); % edit increase minimum area 
conn = 6; % connectivity (i.e. 6-pt)
segment = bwareaopen(segment,areaThresh,conn); % inverse fill holes % note segmentation based on area 

% save raw (cropped) images to imageData structure (for Visual Tool)
imageData.MAG = MAG;
imageData.CD = timeMIP; 
imageData.V = vMean;
imageData.Segmented = segment;

% Feature Extraction - create centerline data 
sortingCriteria = 3; %sorts branches by junctions/intersects 
% spurLength = 15; % 
spurLength = 30; % TEMP: currently very high minimum branch length (removes short spurs) 
[~,~,branchList,~] = feature_extraction(sortingCriteria,spurLength,vMean,segment,handles);

disp('paramMap_params_CSF')
% Return structs instead of numerous variables 
[CLVALS, Planes, Vplanes, CBFSEG, CSFSEG, StructCS, CSFROI, METADATA] ... 
    = paramMap_params_CSF(METADATA,branchList,CD,vMean,MAG,IDX,handles,VOLS); 
METADATA.PCVIPR_HEADER = PCVIPR_HEADER; 
disp('paramMap_params_CSF done')

set(handles.TextUpdate,'String','All Data Loaded'); drawnow;

return





