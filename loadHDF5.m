% function [nframes,matrix,res,timeres,VENC,area_val,diam_val,flowPerHeartCycle_val, ...
%     maxVel_val,PI_val,RI_val,flowPulsatile_val,velMean_val, ...
%     VplanesAllx,VplanesAlly,VplanesAllz,Planes,branchList,segment,r, ...
%     timeMIPcrossection,segmentFull,vTimeFrameave,MAGcrossection, imageData, ...
%     bnumMeanFlow,bnumStdvFlow,StdvFromMean, ... 
%     VplanesCSF, flowCSF, StructCS, CSFSEG, T, CSFROI] = loadHDF5(METADATA.directory,handles,METADATA.infile)

% TODO: fix return
function [CLVALS, Planes, Vplanes, CBFSEG, CSFSEG, StructCS, CSFROI, METADATA, imageData, branchList, flowCSF] = loadHDF5(METADATA, handles)

% LOADHDF5: loadhdf5 reads in PCVIPR data saved in h5
%   Used by: paramMap.m
%   Dependencies: background_phase_correction.m, evaluate_poly.m, calc_cd_timeavg.m,
%   resample_pcvipr_to_isotropic.m (only if METADATA.INTERP), load_vMean.m, feature_extraction.m, slidingThreshold.m
%   METADATA.INTERP: true = isotropic + phase resampling; false = native cropped grid (load_vMean, HDF5 4-D in params)

%% Read HDF5 (update for CBF/CSF combined .h5-files) 
set(handles.TextUpdate,'String','Loading .HDF5 Data'); 
drawnow;

% Resolve input filename robustly (caller may not define METADATA.infile)
flowFile = 'CFLOW.h5';
if isfield(METADATA, 'infile') && ~isempty(METADATA.infile)
    flowFile = METADATA.infile;
end
if ~exist(fullfile(METADATA.directory, flowFile), 'file')
    h5files = dir(fullfile(METADATA.directory, '*.h5'));
    if ~isempty(h5files)
        names = {h5files.name};
        idx = find(strcmpi(names, 'CFLOW.h5'), 1, 'first');
        if isempty(idx)
            idx = 1;
        end
        flowFile = names{idx};
    end
end
METADATA.infile = flowFile;
 
% Read PCVIPR header 
headerFile = [METADATA.directory '/pcvipr_header.txt']; 
PCVIPR_HEADER = struct();
if exist(headerFile, 'file')
    try
        PCVIPR_HEADER = read_pcvipr_header(headerFile);
    catch ME
        warning('loadHDF5:HeaderReadFailed', ...
            'Could not read pcvipr_header.txt (%s). Using default VENC.', ME.message);
    end
else
    warning('loadHDF5:HeaderMissing', ...
        'pcvipr_header.txt not found. Using default VENC.');
end

% Optional header fallback: keep pipeline running without pcvipr_header.txt
% Assumed defaults when header is missing:
%   CSF VENC = 1 cm/s (10 mm/s)
%   CBF VENC = 100 cm/s (1000 mm/s)
if isfield(PCVIPR_HEADER, 'VENC') && ~isempty(PCVIPR_HEADER.VENC)
    METADATA.VENC_CBF_cmps = PCVIPR_HEADER.VENC;
else
    METADATA.VENC_CBF_cmps = 100;
end
METADATA.VENC_CSF_cmps = 1;
METADATA.VENC_CBF_mmps = 10 * METADATA.VENC_CBF_cmps;
METADATA.VENC_CSF_mmps = 10 * METADATA.VENC_CSF_cmps;

% Angiogram is based on blood-flow velocity
vencForAngio = METADATA.VENC_CBF_cmps;

% Frame count fallback when header metadata is unavailable
if isfield(PCVIPR_HEADER, 'nframes') && ~isempty(PCVIPR_HEADER.nframes)
    METADATA.nframes = PCVIPR_HEADER.nframes;
elseif ~isfield(METADATA, 'nframes') || isempty(METADATA.nframes)
    METADATA.nframes = 20;
end

% Standard waveform interpolation resolution (set_metadata also assigns 1000 after load)
METADATA.iframes = 1000;

% add all anatomical information to a volume struct 
StructVols = []; 
StructVols.mcbf = h5read(fullfile(METADATA.directory, flowFile),'/bmag');
StructVols.mcsf = h5read(fullfile(METADATA.directory, flowFile),'/cmag');
StructVols.scsf = h5read(fullfile(METADATA.directory, flowFile),'/bstd');
StructVols.scbf = h5read(fullfile(METADATA.directory, flowFile),'/cstd');
StructVols.cube = ones(size(StructVols.mcsf), 'like', StructVols.mcsf); % fallback when CUBE is unavailable
% StructVols.rage = h5read(fullfile(METADATA.directory, METADATA.infile),'/rage');
% StructVols.pvas = h5read(fullfile(METADATA.directory, METADATA.infile),'/pvasc');
% StructVols.gcsf = h5read(fullfile(METADATA.directory, METADATA.infile),'/gcsf');

% Calculate start and end inds for cropping 3D and 4D data
MAG = single(StructVols.mcbf);
IDX = calc_crop(MAG);

% Crop MAG
MAG = MAG(IDX.start(1):IDX.end(1), IDX.start(2):IDX.end(2), IDX.start(3):IDX.end(3));

% Crop all other 3D volumes (same bounding box)
fns = fields(StructVols);
for j = 1:numel(fns)
    fn = fns{j};
    StructVols.(fn) = StructVols.(fn)(IDX.start(1):IDX.end(1), IDX.start(2):IDX.end(2), IDX.start(3):IDX.end(3));
end

useInterp = isfield(METADATA, 'INTERP') && ~isempty(METADATA.INTERP) && logical(METADATA.INTERP(1));
vel4dResampled = [];

if useInterp
    % Cropped 4-D phase velocity (÷1000, same as load_vMean). CBF = blood (/bph*); CSF = /cph*
    xs = IDX.start(1); xe = IDX.end(1);
    ys = IDX.start(2); ye = IDX.end(2);
    zs = IDX.start(3); ze = IDX.end(3);
    fp = fullfile(METADATA.directory, flowFile);
    rd = @(p) single(h5read(fp, p));
    vel4dNative = struct();
    % CBF (blood)
    vel4dNative.bphx = rd('/bphx') / 1000; vel4dNative.bphx = vel4dNative.bphx(xs:xe, ys:ye, zs:ze, :);
    vel4dNative.bphy = rd('/bphy') / 1000; vel4dNative.bphy = vel4dNative.bphy(xs:xe, ys:ye, zs:ze, :);
    vel4dNative.bphz = rd('/bphz') / 1000; vel4dNative.bphz = vel4dNative.bphz(xs:xe, ys:ye, zs:ze, :);
    % CSF
    vel4dNative.cphx = rd('/cphx') / 1000; vel4dNative.cphx = vel4dNative.cphx(xs:xe, ys:ye, zs:ze, :);
    vel4dNative.cphy = rd('/cphy') / 1000; vel4dNative.cphy = vel4dNative.cphy(xs:xe, ys:ye, zs:ze, :);
    vel4dNative.cphz = rd('/cphz') / 1000; vel4dNative.cphz = vel4dNative.cphz(xs:xe, ys:ye, zs:ze, :);

    if ~isfield(METADATA, 'fov_mm') || isempty(METADATA.fov_mm)
        METADATA.fov_mm = [40, 240, 240];
    end
    if ~isfield(METADATA, 'nframes_target') || isempty(METADATA.nframes_target)
        METADATA.nframes_target = 20;
    end
    METADATA.targetIsotropic_mm = 0.5;

    nfT = METADATA.nframes_target(1);
    tgt = METADATA.targetIsotropic_mm(1);
    set(handles.TextUpdate, 'String', sprintf( ...
        'Resampling to %.2f mm isotropic + %d cardiac frames (imresize3)', tgt, nfT)); drawnow;
    [StructVols, MAG, vMean, METADATA, vel4dResampled] = resample_pcvipr_to_isotropic(StructVols, MAG, vel4dNative, METADATA);
    clear vel4dNative
else
    set(handles.TextUpdate, 'String', 'Native resolution (METADATA.INTERP = false)'); drawnow;
    newDIM = size(MAG);
    vMean = load_vMean(newDIM, METADATA, IDX);
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

% Complex-difference magnitude from time-averaged MAG + vMean (not a max over time).
set(handles.TextUpdate,'String','Computing CD magnitude'); drawnow;
cdImg = calc_cd_timeavg(MAG, vMean, vencForAngio);

% Automatic mask (fallback) + masked CD for interactive review
[segmentFallback, cdMasked] = segment_angiogram_centerline(cdImg, MAG);
set(handles.TextUpdate, 'String', 'Review centerline mask — Export and continue when ready'); drawnow;
% Blocks until Export (chosen mask) or window close (uses automatic mask)
segment = preview_centerline_segmentation(cdMasked, segmentFallback, MAG);

% save raw (cropped) images to imageData structure (for Visual Tool)
imageData.MAG = MAG;
imageData.CD = cdImg;
imageData.V = vMean;
imageData.segment = segment;

% Feature Extraction - create centerline data 
sortingCriteria = 3; %sorts branches by junctions/intersects 
% MinBranchLength for bwskel: lower keeps short intracranial twigs; higher removes spurs/noise
spurLength = 12;
[~,~,branchList,~] = feature_extraction(sortingCriteria,spurLength,vMean,segment,handles);

disp('paramMap_params_CSF')
% Return structs instead of numerous variables 
[CLVALS, Planes, Vplanes, CBFSEG, CSFSEG, StructCS, CSFROI, METADATA, flowCSF] ...
    = paramMap_params_CSF(METADATA, branchList, cdImg, vMean, MAG, IDX, handles, StructVols, vel4dResampled);
METADATA.PCVIPR_HEADER = PCVIPR_HEADER; 
disp('paramMap_params_CSF done')

set(handles.TextUpdate,'String','All Data Loaded'); drawnow;

return





