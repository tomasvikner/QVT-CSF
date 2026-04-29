function varargout = paramMap(varargin)
% TEMP: no dscription 

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @paramMap_OpeningFcn, ...
    'gui_OutputFcn',  @paramMap_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end % End initialization code - DO NOT EDIT

% --- Executes just before paramMap is made visible.
function paramMap_OpeningFcn(hObject, ~, handles, varargin)
disp('paramMap QVT-CSF -- Reminder to load case: ');

% METADATA.INTERP (logical): true = loadHDF5 isotropic + phase resampling; false = native grid / load_vMean path
global METADATA
if isempty(METADATA) || ~isstruct(METADATA)
    METADATA = struct();
end
if ~isfield(METADATA, 'INTERP') || isempty(METADATA.INTERP)
    METADATA.INTERP = false;
end

% add new handles and positions to the GUI object 
axesHandles = findall(hObject, 'Type', 'axes');

handles.CSF1 = axesHandles(4);
handles.CSF2 = axesHandles(1);
handles.CSF3 = axesHandles(2);
handles.CSF4 = axesHandles(3);

% Choose default command line output for paramMap
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
movegui(handles.ParameterTool,'northeast'); %move to top left (WORK)
set(handles.TextUpdate,'String','Load in a 4D Flow Dataset');

% --- Executes on button press in LoadData.
function LoadData_Callback(hObject, ~, handles)
handles.output = hObject; % Choose default command line output for paramMap
guidata(hObject, handles); % Update handles structure

% Create global namespace
global branchList Planes hfull hplane 
global areaThresh fullCData segment
global CLVALS METADATA CBFSEG CSFSEG Vplanes StructCS CSFROI flowCSF
global dcm_obj fig hpatch hscatter cbar hDataTip

% Initial Variables
hfull = handles;
METADATA.Ntxt = []; % used in cursor updatefunction
directory = uigetdir; % interactive directory selection
METADATA.directory = directory; % needed by loadHDF5 when loading from scratch
METADATA.infile = 'CFLOW.h5'; % default combined flow file in selected directory
if ~isfield(METADATA, 'INTERP') || isempty(METADATA.INTERP)
    METADATA.INTERP = false; % see paramMap_OpeningFcn; set true before load for isotropic resampling
end

% Creates list of all .mat files in selected directory
d = dir([directory filesep '*.mat']);
names = {d.name};
isCaseMat = ~endsWith(names, '-ds.mat') & ~endsWith(names, '-vp.mat');
caseNames = names(isCaseMat);

% Also expose preprocessed cases saved as only *-ds.mat/*-vp.mat pairs
isDsMat = endsWith(names, '-ds.mat');
dsNames = names(isDsMat);
for i = 1:numel(dsNames)
    base = dsNames{i}(1:end-7); % strip '-ds.mat'
    caseNames{end+1} = [base '.mat']; %#ok<AGROW>
end

% Keep unique entries while preserving original order
[~, ia] = unique(caseNames, 'stable');
caseNames = caseNames(sort(ia));
fn = [{'Load New Case'}, caseNames];
[fileIndx,~] = listdlg('PromptString','Select a file:', ...
    'ListSize',[200 300],'SelectionMode','single','ListString',fn);
disp(['fileIndx: ' num2str(fileIndx)])

% Data Loading
if fileIndx > 1 % if a pre-processed case is selected
    
    set(handles.TextUpdate,'String','Loading Preprocessed Data'); drawnow;
    METADATA.caseFilePath = [directory filesep fn{fileIndx}];
    
    % Separately load larger than usual data struct 
    ds_filepath = [METADATA.caseFilePath(1:end-4), '-ds.mat'];
    if ~exist(ds_filepath, 'file')
        error('Preprocessed data file missing: %s', ds_filepath);
    end
    load(ds_filepath,'DS')
    vp_filepath = [METADATA.caseFilePath(1:end-4), '-vp.mat'];
    if ~exist(vp_filepath, 'file')
        error('Preprocessed Vplanes file missing: %s', vp_filepath);
    end
    load(vp_filepath, 'Vplanes'); 

    disp('Unpacking values from DS...'); 
    [METADATA, CLVALS, branchList, CBFSEG, Planes, StructCS, CSFSEG, CSFROI, segment, flowCSF] = unpack_struct(DS);
    if isempty(branchList)
        % Backward-compat fallback for older/incomplete preprocessed files
        n = 0;
        if ~isempty(Planes)
            n = size(Planes,1);
            xyz = squeeze(mean(Planes,2)); % center of each plane from 4 corners
            if size(xyz,2) >= 3
                branchList = [xyz(:,1:3), ones(n,1)];
            end
        end
        if isempty(branchList) && ~isempty(CBFSEG)
            n = size(CBFSEG,1);
            branchList = [zeros(n,3), ones(n,1)];
        end
        if isempty(branchList) && isfield(CLVALS, 'flowTA') && ~isempty(CLVALS.flowTA)
            n = numel(CLVALS.flowTA);
            branchList = [zeros(n,3), ones(n,1)];
        end
    end

    set(handles.TextUpdate,'String','Loading Complete'); drawnow; pause(1); 
    set(handles.TextUpdate,'String','Please Select Analysis Plane Location'); drawnow;

else % Load data from scratch
    
    disp('Loading HDF5') % This in turn calls paramMap_params_CSF 
    [CLVALS, Planes, Vplanes, CBFSEG, CSFSEG, StructCS, CSFROI, METADATA, imageData, branchList, flowCSF] = loadHDF5(METADATA, handles); 
    segment = imageData.segment; clear imageData; % maybe not need all "imageData"
    
    time = datestr(now); %#ok<TNOW1,DATST>
    saveState = [time(1:2) time(4:6) time(10:11) '_' time(13:14) time(16:17)];
    set(handles.TextUpdate,'String',['Saving Data as Flow4D_' saveState '.mat']); drawnow;

    disp('Assigning values to DS...'); 
    DS = pack_struct(METADATA, CLVALS, branchList, CBFSEG, Planes, CSFSEG, CSFROI, segment, flowCSF); % Vplanes too large; save alone 
    
    % Saves processed data in same location as pcvipr.mat files
    METADATA.caseFilePath = fullfile(directory,['Flow4D_' saveState '.mat']);

    ds_filepath = [METADATA.caseFilePath(1:end-4), '-ds.mat'];
    save(ds_filepath, 'DS', '-v7.3') % if too large 
    vp_filepath = [METADATA.caseFilePath(1:end-4), '-vp.mat'];
    save(vp_filepath, 'Vplanes', '-v7.3')

end
METADATA = set_metadata(METADATA); % additional settings not loaded
if ~isfield(METADATA, 'nframes') || isempty(METADATA.nframes)
    METADATA.nframes = 20;
end
if ~isfield(METADATA, 'coregTarget') || isempty(METADATA.coregTarget)
    METADATA.coregTarget = 'cnob';
end
METADATA.vesselsAnalyzed = {};

% Makes directory if it does already exist (folder is time-stamped)
warning off
METADATA.SavePath = [directory '/OUTPUT/']; % edit 
if ~exist(METADATA.SavePath) %#ok<*EXIST>
    mkdir(METADATA.SavePath)
    disp(['mkdir: ' METADATA.SavePath]); 
end

%%% Plotting 3D Interactive Display
set(handles.parameter_choice,'Value',3); %set parameter to flow as default
set(handles.Transparent, 'Value',0);
set(handles.AreaThreshSlide, 'Value',0);

% Initialize visualization
fig = figure(1); cla
isoAx = [];
if isfield(hfull, 'isoAx') && ~isempty(hfull.isoAx) && isgraphics(hfull.isoAx)
    isoAx = hfull.isoAx;
elseif isfield(handles, 'isoAx') && ~isempty(handles.isoAx) && isgraphics(handles.isoAx)
    isoAx = handles.isoAx;
else
    isoAx = gca;
end
if isempty(segment)
    set(handles.TextUpdate,'String','Loaded preprocessed data (no 3D segment found)'); drawnow;
else
    hpatch = show_iso_seg(isoAx, segment, hpatch);
end

% Dark background for 3D centerline figure (defaults are often white on newer MATLAB)
if isgraphics(isoAx)
    axes(isoAx);
    set(isoAx, 'Color', [0 0 0], 'XColor', [0.75 0.75 0.75], 'YColor', [0.75 0.75 0.75], 'ZColor', [0.75 0.75 0.75]);
end
set(fig, 'Color', [0 0 0]);

% Turn on data cursormode within the figure
dcm_obj = datacursormode(fig); %create dataCursorManager object
datacursormode on;
dcm_obj.DisplayStyle = 'window';
set(handles.CBARmin,'String','min')
set(handles.CBARmax,'String','max')
METADATA.branchLabeled = 0;

% This will be used in the update function for cursor text
METADATA.Labeltxt = {'Flow: ',  ' mL/s ';'Average: ',' mL/s '};
cdata = CLVALS.flowTA;

% Normalize branchList shape for robust reload of preprocessed cases
if isvector(branchList)
    if mod(numel(branchList), 4) == 0
        branchList = reshape(branchList, [], 4);
    elseif mod(numel(branchList), 3) == 0
        branchList = reshape(branchList, [], 3);
    end
end
if size(branchList,2) < 3 && size(branchList,1) >= 3
    branchList = branchList';
end
if size(branchList,2) < 3
    error('Loaded branchList has invalid size [%d %d]; expected at least 3 columns.', ...
        size(branchList,1), size(branchList,2));
end
if isrow(cdata)
    cdata = cdata';
end
if numel(cdata) ~= size(branchList,1)
    cdata = cdata(:);
    n = min(numel(cdata), size(branchList,1));
    branchList = branchList(1:n,:);
    cdata = cdata(1:n);
end
hold on
hscatter = scatter3(branchList(:,1),branchList(:,2),branchList(:,3),METADATA.dotSize,cdata,'filled');
hold off

clim([min(cdata) max(cdata)]);
cbar = colorbar;
clim([0 0.8*max(CLVALS.flowTA(:))])
set(get(cbar,'xlabel'),'string','Flow (mL/s)','fontsize',16,'Color','white');
set(cbar,'FontSize',16,'color','white');
ax = gca;
r = 10;
if isfield(METADATA, 'PLANESIZE') && ~isempty(METADATA.PLANESIZE)
    r = METADATA.PLANESIZE;
end
xlim([ax.XLim(1)-r ax.XLim(2)+r]) %buffer with extra space for planes
ylim([ax.YLim(1)-r ax.YLim(2)+r])
zlim([ax.ZLim(1)-r ax.ZLim(2)+r])

% Initialize visualization of tangent planes
hold on
hplane = fill3(Planes(1,:,1)',Planes(1,:,2)',Planes(1,:,3)',[1 0 0], ...
    'EdgeColor',[1 0 0],'FaceAlpha',0.3,'PickableParts','none', ...
    'Parent', fig.CurrentAxes); % fill3(pty',ptx',ptz','r') for isosurface
hold off

% Update string (undocumentedmatlab.com/articles/controlling-plot-data-tips)
set(dcm_obj,'UpdateFcn',@myupdatefcn_all); % update dataCursor w/ cust. fcn
hDataTip = dcm_obj.createDatatip(hscatter); % NOTE: is this used? 

% Convert toolbar to old style and add hot keys
fig.CurrentAxes.Toolbar = [];
addToolbarExplorationButtons(fig)
updateDataCursors(dcm_obj)

% Calculate average area per branch
areaThresh.AvgArea = size(max(branchList(:,4)),1);
for n = 1:max(branchList(:,4))
    Btemp = branchList(:,4)==n;
    areaThresh.AvgArea(n,1) = mean(CLVALS.area(Btemp)); %mean area of branch
end

areaThresh.LogPoints = true(size(branchList,1),1); % logical array of 1s for areaThresh
fullCData = CLVALS.flowTA; % initialize fullCData color as flow

nf = max(2, METADATA.nframes);
steps = [1./(nf-1) 10./(nf-1)]; % set so one 'slide' moves to the next slice exactly
set(handles.VcrossTRslider,'SliderStep',steps);

% new initializations to track manual coregs and segmentations 
nplanes = size(branchList, 1);
CSFSEG.coregTrack = zeros(nplanes, 1);
CSFSEG.mansegTrack = zeros(nplanes, 1);
METADATA.wfState = 'PC-1'; % start at PCA waveform
disp(['LoadData_Callback complete; WF state: ' METADATA.wfState])

% --- Outputs from this function are returned to the command line.
function varargout = paramMap_OutputFcn(~, ~, handles)
varargout{1} = handles.output;

% --- Executes on selection change in parameter_choice.
function parameter_choice_Callback(~, ~, handles)
global hscatter METADATA fig cbar dcm_obj fullCData 

disp('Called parameter_choice_Callback')

% Parameter definitions from CLVALS 
P = get_parameter_table(CLVALS);

str   = get(handles.parameter_choice,'String');
param = str{get(handles.parameter_choice,'Value')};
D     = P.(param);

CData    = D{1};
clim     = D{2};
METADATA.Labeltxt = D{3};
cb_label = D{4};

% Auto clim if not specified
if isempty(clim), clim = [min(CData) max(CData)]; end
hscatter.CData = CData;
fullCData      = CData;
clim(fig.CurrentAxes, clim)
set(get(cbar,'xlabel'),'String',cb_label,'FontSize',16,'Color','white')
set(cbar,'FontSize',16,'Color','white')
set(handles.CBARmin,'String',num2str(clim(1)),'Value',clim(1))
set(handles.CBARmax,'String',num2str(clim(2)),'Value',clim(2))
updateDataCursors(dcm_obj)

% --- Executes during object creation, after setting all properties.
function parameter_choice_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function plot_flowWaveform_CreateFcn(~, ~, ~) % note this function is never called and does nothing? 
disp('Called plot_flowWaveform_CreateFcn; doing nothing?')

% --- Executes on slider movement.
function Transparent_Callback(hObject, ~, ~)
global hpatch Sval
Sval = get(hObject,'Value');
set(hpatch,'FaceAlpha',Sval);

% --- Executes during object creation, after setting all properties.
function Transparent_CreateFcn(hObject, ~, ~)
global Sval
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
set(hObject, 'Min', 0);
set(hObject, 'Max', 1);
set(hObject,'Value',0);
Sval = get(hObject,'Value');

function CBARmin_Callback(~, ~, handles)
global fig
maxV = str2double(get(handles.CBARmax,'String'));
minV = str2double(get(handles.CBARmin,'String'));
clim(fig.CurrentAxes,[minV maxV])

% --- Executes during object creation, after setting all properties.
function CBARmin_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function CBARmax_Callback(~, ~, handles)
global fig
maxV =   str2double(get(handles.CBARmax,'String'));
minV =   str2double(get(handles.CBARmin,'String'));
clim(fig.CurrentAxes,[minV maxV])

% --- Executes during object creation, after setting all properties.
function CBARmax_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in CBARselection.
function CBARselection_Callback(hObject, ~, ~)
global fig
contents = cellstr(get(hObject,'String')); %turn color options to cells
colormap(fig.Children(end),contents{get(hObject,'Value')})

% --- Executes during object creation, after setting all properties.
function CBARselection_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in NamePoint.
function NamePoint_Callback(~, ~, handles)
global METADATA dcm_obj 

contents = cellstr(get(handles.NamePoint,'String'));
METADATA.PointLabel = contents{get(handles.NamePoint,'Value')};
if sum(contains(METADATA.vesselsAnalyzed,METADATA.PointLabel))
    set(handles.NamePoint,'ForegroundColor',[0.6 0.6 0.6]);
else
    set(handles.NamePoint,'ForegroundColor',[0 0 0]);
end 
updateDataCursors(dcm_obj)

% --- Executes during object creation, after setting all properties.
function NamePoint_CreateFcn(hObject, ~, ~)
global METADATA
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
contents = cellstr(get(hObject,'String'));
METADATA.PointLabel = contents{get(hObject,'Value')};

% --- Executes while toggling waveform/segmentation status. 
function toggleCnoB_Callback(~, ~, ~)
global METADATA dcm_obj
if ~strcmp(METADATA.wfState, 'CnoB')
    disp('WF switch to CnoB');
    METADATA.wfState = 'CnoB';
end
set(dcm_obj,'UpdateFcn',@myupdatefcn_all); 

function togglePCA_Callback(~, ~, ~)
global METADATA dcm_obj
if ~strcmp(METADATA.wfState, 'PC-1')
    disp('WF switch to PC-1');
    METADATA.wfState = 'PC-1';
end
set(dcm_obj,'UpdateFcn',@myupdatefcn_all); 

function toggleCUBE_Callback(~, ~, ~)
global METADATA dcm_obj
if ~strcmp(METADATA.wfState, 'CUBE')
    disp('WF switch to CUBE');
    METADATA.wfState = 'CUBE';
end
set(dcm_obj,'UpdateFcn',@myupdatefcn_all); 

% --- Executes on button press: 
function manualSeg_Callback(~, ~, ~)
global METADATA flowCSF dcm_obj CSFSEG CBFSEG branchList 
mask = calc_manseg(dcm_obj, branchList, CBFSEG); 
[pindex, ~] = get_pindex(dcm_obj, branchList);
CSFSEG.(METADATA.coregTarget)(pindex,:) = mask(:); % Update segmentation globally 
CSFSEG.segTrack(pindex) = 1;
flowCSF = updateWaveforms(flowCSF, METADATA.coregTarget, VplanesCSF, CSFSEG, ...
    pindex, METADATA.PLANSIZE, METADATA.nframes); % Update CSF WFs 
set(dcm_obj, 'UpdateFcn', @myupdatefcn_all);

% --- Executes on button press in SavePoint.
function SavePoint_Callback(handles, ~, ~)
global CLVALS METADATA CSFROI branchList dcm_obj WFPS VcrossTR CcrossTR CC 
METADATA.vesselsAnalyzed{end+1} = METADATA.PointLabel;
set(handles.TextUpdate,'String','Saving Data.'); drawnow;
info_struct = getCursorInfo(dcm_obj);
save_point(info_struct, branchList, CLVALS, VcrossTR, CcrossTR, CC, CSFROI, WFPS, METADATA.PointLabel, METADATA) % Cleaner save function 

% --- NoteBox_Callback
function NoteBox_Callback(~, ~, ~)

% --- Executes during object creation, after setting all properties.
function NoteBox_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in AxialView.
function AxialView_Callback(~, ~, ~)
global fig; view(fig.CurrentAxes,[180,90])

% --- Executes on button press in SagittalView.
function SagittalView_Callback(~, ~, ~)
global fig; view(fig.CurrentAxes,[180,0])

% --- Executes on button press in CoronalView.
function CoronalView_Callback(~, ~, ~)
global fig; view(fig.CurrentAxes,[90,0])

% --- Executes on slider movement.
function AreaThreshSlide_Callback(hObject, handles, ~)
global areaThresh branchList hscatter CLVALS
areaThresh.LogPoints = find(areaThresh.AvgArea>max(areaThresh.AvgArea)*get(hObject,'Value')*.15);
areaThresh.LogPoints = ismember(branchList(:,4),areaThresh.LogPoints);
hscatter = area_thresh(handles, hscatter, CLVALS, branchList, areaThresh.LogPoints);

% --- Executes during object creation, after setting all properties.
function AreaThreshSlide_CreateFcn(hObject, ~, ~)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% --- Executes on slider movement.
function clWidthSlider_Callback(hObject, ~, ~)
global hscatter; set(hscatter,'SizeData',get(hObject,'Value'));

% --- Executes during object creation, after setting all properties.
function clWidthSlider_CreateFcn(hObject, ~, ~)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% --- Executes on button press in VisualTool.
function VisualTool_Callback(handles, ~, ~)
global Planes branchList segment METADATA res
set(handles.TextUpdate,'String','Opening Visual Tool'); drawnow;
fourDvis(Planes,branchList,segment,METADATA);
uiwait;
set(handles.TextUpdate,'String','Visual Tool Closed'); drawnow;

% --- Executes on button press in InvertArea.
function InvertArea_Callback(hObject, handles, ~)
global branchList hscatter CLVALS
hscatter = invert_area(hObject, hscatter, branchList, CLVALS, handles); 

% --- Executes on slider movement.
function VcrossTRslider_Callback(handles, ~, ~)
updateVcrossTR(handles)

% --- Executes during object creation, after setting all properties.
function VcrossTRslider_CreateFcn(hObject, ~, ~)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function updateVcrossTR(~, ~, ~)
global dcm_obj hfull CBFSEG Vplanes CC CSFSEG %#ok<*GVMIS> % add some TR CSF planes 
global METADATA branchList 

info_struct = getCursorInfo(dcm_obj);
if isempty(info_struct)
    return;
end

    [pindex, ~] = get_pindex(dcm_obj, branchList);
    pindex = pindex(1);
    imdim = sqrt(size(CBFSEG,2)); % side length of cross-section
    Maskcross = reshape(CBFSEG(pindex,:), imdim, imdim);

    % Get time-resolved blood (V) and CSF (C) 
    [VcrossTR, CcrossTR] = get_vCS_TR(Vplanes, pindex, imdim, METADATA);
    nTr = size(VcrossTR, 3);
    sliceNum = 1 + round(get(hfull.VcrossTRslider,'Value') * max(0, nTr - 1));
    sliceNum = max(1, min(sliceNum, nTr));
    Vsl = VcrossTR(:, :, sliceNum);
    [CC, ~] = reshape_all_CS(CSFSEG, CBFSEG, pindex);

    show_img_and_roi(hfull.TRcross, Vsl, ...
        'CBF velocity (time-res.)', Maskcross, [min(Maskcross(:).*Vsl(:)), max(Maskcross(:).*Vsl(:))]);

    show_img_and_roi(hfull.CSF4, CcrossTR(:,:,sliceNum), ...
        'CSF velocity (time-res.)', CC.madj, [min(CcrossTR(:)) max(CcrossTR(:))]);

    if strcmp(METADATA.wfState, 'CUBE')
        visboundaries(hfull.CSF4,CC.bcube,'LineWidth',1) % auto/manual combo to vis/seg
    else
        visboundaries(hfull.CSF4,CC.madj,'LineWidth',1) % auto/manual combo to vis/seg
    end

function txt = myupdatefcn_all(~, ~, ~) 

% Avoid close(2:3): it can destroy unrelated figures; manual seg closes its own fig.

% Customizes text of data tips
% TODO: fullCData, METADATA.branchLabeled? 
global METADATA branchList fullCData Planes hplane dcm_obj hfull flowCSF CBFSEG Vplanes CSFSEG CC WFPS CLVALS fig

fprintf(1, '[myupdatefcn_all] datacursor update\n');

ax = hfull.pfwaveform;
txt = {''};

try
% --- Extract inds at pointer (pointer, range, range-rearranged) 
[pindex, index_range, pinds] = get_pindex(dcm_obj, branchList); 
pindex = pindex(1);

% --- Update Cross-sectional planes for points
set(hplane,'XData',Planes(pindex,:,1)','YData',Planes(pindex,:,2)','ZData',Planes(pindex,:,3)')

% --- Reshape local CS for vis 
[CC, imdim] = reshape_all_CS(CSFSEG, CBFSEG, pindex); 

% --- Extract time-resolved CBF and CSF velocities from Vplanes and pindex 
[VcrossTR, CcrossTR] = get_vCS_TR(Vplanes, pindex, imdim, METADATA);
nTr = size(VcrossTR, 3);
sliceNum = 1 + round(get(hfull.VcrossTRslider,'Value') * max(0, nTr - 1));
sliceNum = max(1, min(sliceNum, nTr));
Vsl = VcrossTR(:, :, sliceNum);

% --- Visualize CBF 
Maskcross = reshape(CBFSEG(pindex,:), imdim, imdim);
show_img_and_roi(hfull.MAGcross, CC.mcbf, 'CBF magnitude', Maskcross);
if isfield(CC, 'rage')
    show_img_and_roi(hfull.CDcross, CC.rage,'T1-w MP-RAGE', Maskcross);
else
    show_img_and_roi(hfull.CDcross, CC.CD,'Complex Difference', Maskcross);
end
show_img_and_roi(hfull.VELcross, CC.CD, 'Complex Difference', Maskcross);
show_img_and_roi(hfull.TRcross, Vsl, ...
    'CBF velocity (time-res.)', Maskcross, [min(Maskcross(:).*Vsl(:)), max(Maskcross(:).*Vsl(:))]);

% --- Visualize CSF 
show_img_and_roi(hfull.CSF1, CC.mcsf, 'CSF magnitude', CC.bics);
show_img_and_roi(hfull.CSF2, CC.cube, 'CUBE Anti-FLAIR', CC.bcube);
show_img_and_roi(hfull.CSF3, CC.scsf, 'CSF velocity STD', CC.bcsf );
show_img_and_roi(hfull.CSF4, CcrossTR(:,:,sliceNum), ...
'CSF velocity (time-res.)', CC.madj, [min(CcrossTR(:)) max(CcrossTR(:))]);

% --- Calculate median WFs in center points and WFs in surrounding points 
[CBF, CSF, cbfwf, csfwf] = calc_waveforms(CLVALS, flowCSF, pinds, METADATA.wfState);

% --- Interp prep before coupling analysis
[CBF, CSF, cbfwf, csfwf] = interpCoupling(METADATA, CBF, CSF, cbfwf, csfwf);
cardiacCycle = 1:numel(CBF);

% --- Calculate WF (CBF-CSF) coupling
[rmax, mlag] = waveformCoupling(CBF, CSF, METADATA.MAXLAG); % max 300 ms delay?
ttext = ['Coupling (xcorr [-1, 1] / lag (ms): ' num2str(rmax, 3), ' / ' num2str(mlag, 3)];

% --- Amplitude and stroke volumes 
[amp, cdv] = calc_metrics(CSF, CBF);

% --- Pack into waveform point save struct (explicit fields: pack_struct needs named
%     variables; CBF' / CSF' are expressions and break inputname inside pack_struct)
WFPS = struct();
WFPS.amp = amp;
WFPS.cdv = cdv;
WFPS.rmax = rmax;
WFPS.mlag = mlag;
WFPS.CBF = struct('all', cbfwf, 'avg', CBF(:).');
WFPS.CSF = struct('all', csfwf, 'avg', CSF(:).');
WFPS.irange = index_range;
WFPS.pindex = pindex;

% --- Clear CBF and CSF axes 
yyaxis(hfull.pfwaveform, 'left');
cla(hfull.pfwaveform);
yyaxis(hfull.pfwaveform, 'right');
cla(hfull.pfwaveform);

% --- Plot CSF waveform 
yyaxis(hfull.pfwaveform, 'left');
if METADATA.PATCH
    [X, Y] = calc_patch(WFPS.CSF.all, cardiacCycle);
    plot_waveform_patched(hfull.pfwaveform, X, Y, CSF, cardiacCycle, METADATA.c1, METADATA.PLW);
else
    plot_waveform_lines(ax, WFPS.CSF, cardiacCycle, METADATA.c1, METADATA.PLW)
end
ylabel(ax, 'Flow (CSF) (mL/s)', 'FontSize', 16);

% --- Plot CBF waveform 
yyaxis(hfull.pfwaveform, 'right');
if METADATA.PATCH
    [X, Y] = calc_patch(WFPS.CBF.all, cardiacCycle);
    plot_waveform_patched(hfull.pfwaveform, X, Y, CBF, cardiacCycle, METADATA.c2, METADATA.PLW);
else
    plot_waveform_lines(ax, WFPS.CBF, cardiacCycle, METADATA.c2, METADATA.PLW)
end
ylabel(ax, 'Flow (CBF) (mL/s)', 'FontSize', 16);

% --- Plot CBF/CSF legend and axis 
legend(hfull.pfwaveform, '', ...
['CSF amp.: ' num2str(amp.CSF, 3) ' mL/s'], '', ...
['CBF amp.: ' num2str(amp.CBF, 3) ' mL/s'], ...
'Box', 'off', 'FontSize', 16, 'FontWeight', 'bold', ...
'Location', 'North');
title(hfull.pfwaveform, ttext, 'FontSize', 16);
grid(hfull.pfwaveform, 'on');
hfull.pfwaveform.XAxisLocation = 'bottom';
yyaxis(hfull.pfwaveform, 'left');
axis(hfull.pfwaveform, 'tight');

% --- Labels for centerline visualization
bnum = branchList(pindex,4);
value = fullCData(pindex);  
average = fullCData(index_range);
[METADATA, txt] = calc_labels(METADATA, bnum, branchList, value, average, pindex);

catch ME
    fprintf(2, '[myupdatefcn_all] ERROR: %s\n', ME.message);
    txt = {ME.message};
end

% --- Executes when user attempts to close ParameterTool.
function ParameterTool_CloseRequestFcn(hObject, ~, ~)
delete(hObject);

