function varargout = paramMap_clean(varargin)
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

% add new handles and positions to the GUI object 
axesHandles = findall(hObject, 'Type', 'axes');

% TODO: why this order 
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
global branchList Planes hfull p branchLabeled Ntxt
global AveAreaBranch LogPoints fullCData segment
global CLVALS METADATA CBFSEG CSFSEG Vplanes StructCS CSFROI 
global dcm_obj fig hpatch hscatter Labeltxt cbar hDataTip
global imageData 

% Initial Variables
hfull = handles;
Ntxt = []; % used in cursor updatefunction
directory = uigetdir; % interactive directory selection

% Creates list of all .mat files in selected directory
d = dir([directory filesep '*.mat']);
fn = [{'Load New Case'},{d.name}];
[fileIndx,~] = listdlg('PromptString','Select a file:', ...
    'ListSize',[200 300],'SelectionMode','single','ListString',fn);
disp(['fileIndx: ' num2str(fileIndx)])

% Data Loading
if fileIndx > 1 % if a pre-processed case is selected
    
    set(handles.TextUpdate,'String','Loading Preprocessed Data'); drawnow;
    METADATA.caseFilePath = [directory filesep fn{fileIndx}];
    
    % Separately load larger than usual data struct 
    ds_filepath = [METADATA.caseFilePath(1:end-4), '-ds.mat'];
    load(ds_filepath,'data_struct') 
    vp_filepath = [METADATA.caseFilePath(1:end-4), '-vp.mat'];
    load(vp_filepath, 'Vplanes'); 

    disp('Unpacking values from data_struct...'); 
    [METADATA, CLVALS, branchList, CBFSEG, Planes, StructCS, CSFSEG, CSFROI] = unpack_struct(data_struct);

    set(handles.TextUpdate,'String','Loading Complete'); drawnow; pause(1); 
    set(handles.TextUpdate,'String','Please Select Analysis Plane Location'); drawnow;

else % Load data from scratch
    
    disp('Loading HDF5') % This in turn calls paramMap_params_CSF where all CL c
    [CLVALS, Planes, Vplanes, CBFSEG, CSFSEG, StructCS, CSFROI, METADATA, imageData] = loadHDF5(METADATA, handles); 
    
    time = datestr(now); %#ok<TNOW1,DATST>
    saveState = [time(1:2) time(4:6) time(10:11) '_' time(13:14) time(16:17)];
    set(handles.TextUpdate,'String',['Saving Data as Flow4D_' saveState '.mat']); drawnow;

    disp('Assigning values to data_struct...'); 
    data_struct = pack_struct(METADATA, CLVALS, branchList, CBFSEG, Planes, StructCS, CSFSEG, CSFROI); % Vplanes too large; save alone 
    
    % Saves processed data in same location as pcvipr.mat files
    METADATA.caseFilePath = fullfile(directory,['Flow4D_' saveState '.mat']);

    ds_filepath = [METADATA.caseFilePath(1:end-4), '-ds.mat'];
    save(ds_filepath, 'data_struct', '-v7.3') % if too large 
    vp_filepath = [METADATA.caseFilePath(1:end-4), '-vp.mat'];
    save(vp_filepath, 'Vplanes', '-v7.3')

end
METADATA.directory = directory;
METADATA = set_metadata(METADATA); % additional settings not loaded
METADATA.vesselAnalyzed = {};

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
hpatch = show_iso_segment(hfull.isoAx, segment, hpatch); % TODO: hfull. WHAT AXIS?

% Turn on data cursormode within the figure
dcm_obj = datacursormode(fig); %create dataCursorManager object
datacursormode on;
dcm_obj.DisplayStyle = 'window';
set(handles.CBARmin,'String','min')
set(handles.CBARmax,'String','max')
branchLabeled = 0;

% This will be used in the update function for cursor text
Labeltxt = {'Flow: ',  ' mL/s ';'Average: ',' mL/s '};
cdata = CLVALS.flowTA;
hold on
hscatter = scatter3(branchList(:,1),branchList(:,2),branchList(:,3),METADATA.dotSize,cdata,'filled');
hold off

clim([min(cdata) max(cdata)]);
cbar = colorbar;
clim([0 0.8*max(CLVALS.flowTA(:))])
set(get(cbar,'xlabel'),'string','Flow (mL/s)','fontsize',16,'Color','white');
set(cbar,'FontSize',16,'color','white');
ax = gca;
xlim([ax.XLim(1)-r ax.XLim(2)+r]) %buffer with extra space for planes
ylim([ax.YLim(1)-r ax.YLim(2)+r])
zlim([ax.ZLim(1)-r ax.ZLim(2)+r])

% Initialize visualization of tangent planes
hold on
p = fill3(Planes(1,:,1)',Planes(1,:,2)',Planes(1,:,3)',[1 0 0], ...
    'EdgeColor',[1 0 0],'FaceAlpha',0.3,'PickableParts','none', ...
    'Parent', fig.CurrentAxes); % fill3(pty',ptx',ptz','r') for isosurface
hold off

% Update string (undocumentedmatlab.com/articles/controlling-plot-data-tips)
set(dcm_obj,'UpdateFcn',@myupdatefcn_all); % update dataCursor w/ cust. fcn
hDataTip = dcm_obj.createDatatip(hscatter);

% Convert toolbar to old style and add hot keys
fig.CurrentAxes.Toolbar = [];
addToolbarExplorationButtons(fig)
updateDataCursors(dcm_obj)

% Calculate average area per branch
AveAreaBranch = size(max(branchList(:,4)),1);
for n = 1:max(branchList(:,4))
    Btemp = branchList(:,4)==n;
    AveAreaBranch(n,1) = mean(CLVALS.area(Btemp)); %mean area of branch
end

LogPoints = true(size(branchList,1),1); % logical array of 1s for areaThresh
fullCData = CLVALS.flowTA; %i nitialize fullCData color as flow

steps = [1./(METADATA.nframes-1) 10./(METADATA.nframes-1)]; % set so one 'slide' moves to the next slice exactly
set(handles.VcrossTRslider,'SliderStep',steps);

% new initializations to track manual coregis and segmentations 
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
global hscatter fig cbar dcm_obj fullCData Labeltxt

disp('Called parameter_choice_Callback')

% Parameter definitions from CLVALS 
P = get_parameter_table(CLVALS);

str   = get(handles.parameter_choice,'String');
param = str{get(handles.parameter_choice,'Value')};
D     = P.(param);

CData    = D{1};
clim     = D{2};
Labeltxt = D{3};
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
global PointLabel dcm_obj METADATA.vesselsAnalyzed

contents = cellstr(get(handles.NamePoint,'String'));
PointLabel = contents{get(handles.NamePoint,'Value')};
if sum(contains(METADATA.vesselsAnalyzed,PointLabel))
    set(handles.NamePoint,'ForegroundColor',[0.6 0.6 0.6]);
else
    set(handles.NamePoint,'ForegroundColor',[0 0 0]);
end 
updateDataCursors(dcm_obj)

% --- Executes during object creation, after setting all properties.
function NamePoint_CreateFcn(hObject, ~, ~)
global PointLabel
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
contents = cellstr(get(hObject,'String'));
PointLabel = contents{get(hObject,'Value')};

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
global CLVALS METADATA CSFROI
global PointLabel branchList dcm_obj METADATA.vesselsAnalyzed 
global WFPS VcrossTR CcrossTR CC 
METADATA.vesselsAnalyzed{end+1} = PointLabel;
set(handles.TextUpdate,'String','Saving Data.'); drawnow;
info_struct = getCursorInfo(dcm_obj);
save_point(info_struct, branchList, CLVALS, VcrossTR, CcrossTR, CC, CSFROI, WFPS, PointLabel, METADATA) % Cleaner save function 

% --- NoteBox_Callback
function NoteBox_Callback(~, ~, ~)

% --- Executes during object creation, after setting all properties.
function NoteBox_CreateFcn(hObject)
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
global LogPoints branchList hscatter AveAreaBranch CLVALS
LogPoints = find(AveAreaBranch>max(AveAreaBranch)*get(hObject,'Value')*.15);
LogPoints = ismember(branchList(:,4),LogPoints);
hscatter = area_thresh(handles, hscatter, CLVALS, branchList, LogPoints);

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
global Planes branchList segment METADATA.caseFilePath res 
set(handles.TextUpdate,'String','Opening Visual Tool'); drawnow;
fourDvis(Planes,branchList,segment,METADATA.caseFilePath,res);
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
global dcm_obj hfull CBFSEG Vplanes CC %#ok<*GVMIS> % add some TR CSF planes 
global METADATA branchList 

info_struct = getCursorInfo(dcm_obj);
if ~isempty(info_struct)

    [pindex, ~] = get_pindex(dcm_obj, branchList);
    imdim = sqrt(size(CBFSEG,2)); % side length of cross-section
    Maskcross = reshape(CBFSEG(pindex,:), imdim, imdim);

    % Get slice number from slider
    sliceNum = 1+round( get(hfull.VcrossTRslider,'Value').*(METADATA.nframes-1) ); 

    % Get time-resolved blood (V) and CSF (C) 
    [VcrossTR, CcrossTR] = get_vCS_TR(Vplanes, pindex);

    minn = min(Maskcross.*VcrossTR,[],'all');
    maxx = max(Maskcross.*VcrossTR,[],'all');
    imshow(VcrossTR(:,:,sliceNum),[minn maxx],'InitialMagnification','fit','Parent',hfull.TRcross)
    title(hfull.TRcross, 'CBF velocity (time-res. )', 'FontSize', 13);
    visboundaries(hfull.TRcross,Maskcross,'LineWidth',1)

    % CSF CS data
    minn = min(CcrossTR(:)); maxx = max(CcrossTR(:));
    imshow(CcrossTR(:,:,sliceNum),[minn maxx],'InitialMagnification','fit','Parent',hfull.CSF4)
    title(hfull.CSF4, 'CSF velocity (time-res. )', 'FontSize', 13);

    if strcmp(METADATA.wfState, 'CUBE')
        visboundaries(hfull.CSF4,CC.bcube,'LineWidth',1) % auto/manual combo to vis/seg
    else
        visboundaries(hfull.CSF4,CC.madj,'LineWidth',1) % auto/manual combo to vis/seg
    end
end 

function txt = myupdatefcn_all(~, ~, ~) 

try 
    close(2:3) % closing manual ROI windows when changing pindex  
catch
end

% Customizes text of data tips
% TODO: MAGcrossection and timeMIPcrossection - replace
global Labeltxt branchLabeled branchList fullCData
global METADATA Planes p dcm_obj Ntxt hfull timeMIPcrossection flowCSF
global CBFSEG MAGcrossection
global Vplanes CSFSEG CC % 
global WFPS % waveform point save struct

% Again extracting p_index
[pindex, index_range] = get_pindex(dcm_obj, branchList); 

% Update Cross-sectional planes for points
set(p,'XData',Planes(pindex,:,1)','YData',Planes(pindex,:,2)','ZData',Planes(pindex,:,3)')
imdim = sqrt(size(CBFSEG,2)); %side length of cross-section

% difining new local CS for vis 
CC = [];
fns = fieldnames(CSFSEG);
for i = 1:numel(fns) 
    fn = fns{i};
    if ~contains(fn, 'Track') % Don't update segmentation and coreg track
        CC.(fn) = reshape(double(CSFSEG.(fn)(pindex,:)),imdim,imdim);
    end
end
if sum(CC.madj(:)==0), CC.madj = CC.auto; end

% Extract time-resolved CBF and CSF velocities from Vplanes and pindex 
[VcrossTR, CcrossTR] = get_vCS_TR(Vplanes, pindex);
sliceNum = 1 + round(get(hfull.VcrossTRslider,'Value') * (METADATA.nframes-1));

% ===================== CBF Visual =====================
Maskcross = reshape(CBFSEG(pindex,:), imdim, imdim);
MAGcross = reshape(MAGcrossection(pindex,:), imdim, imdim);
CDcross = reshape(timeMIPcrossection(pindex,:), imdim, imdim);
show_img_and_roi(hfull.MAGcross, MAGcross, 'CBF magnitude', Maskcross);
show_img_and_roi(hfull.CDcross, CC.rage,'T1-w MP-RAGE', Maskcross);
show_img_and_roi( hfull.VELcross, CDcross, 'Complex Difference', Maskcross);
minn = min(Maskcross .* VcrossTR, [], 'all');
maxx = max(Maskcross .* VcrossTR, [], 'all');
show_img_and_roi(hfull.TRcross, VcrossTR(:,:,sliceNum), ...
    'CBF velocity (time-res.)', Maskcross, [minn maxx] );

% ===================== CSF Visual =====================
show_img_and_roi(hfull.CSF1, CC.mcsf, 'CSF magnitude', CC.bics);
show_img_and_roi(hfull.CSF2, CC.cube, 'CUBE Anti-FLAIR', CC.bcube);
show_img_and_roi(hfull.CSF3, CC.scsf, 'CSF velocity STD', CC.bcsf );
minn = min(CcrossTR(:));
maxx = max(CcrossTR(:));
show_img_and_roi(hfull.CSF4, CcrossTR(:,:,sliceNum), ...
'CSF velocity (time-res.)', CC.madj, [minn maxx] );

% Get value of parameter at point and mean within 5pt window
value = fullCData(pindex); % value/average sent to calc_labels 
average = fullCData(index_range);
pside = index_range;
pside(pside==pindex) = [];
if size(pside, 2) > size(pside, 1), pside = pside'; end
pinds = [pindex; pside]; % NOTE: why this pinds definition? 

% Calculate median WFs in center points and WFs in surrounding points 
[CBF, CSF, cbfwf, csfwf] = calc_waveforms(CLVALS, flowCSF, pinds, METADATA.wfState);

% Interp prep before coupling analysis
METADATA.iframes = 1000;
[CBF, CSF, cbfwf, csfwf] = interpCoupling(METADATA, CBF, CSF, cbfwf, csfwf);
cardiacCycle = 1:numel(CBF);

% Calculate WF (CBF-CSF) coupling and visualize
[rmax, mlag] = waveformCoupling(CBF, CSF, METADATA.MAXLAG); % max 300 ms delay?
ttext = ['Coupling (xcorr [-1, 1] / lag (ms): ' num2str(rmax, 3), ' / ' num2str(mlag, 3)];

% Amplitude and stroke volumes 
[amp, cdv] = calc_metrics(CSF, CBF);

% Pack into waveform point save struct 
WFPS = pack_struct(amp, cdv, rmax, mlag, cbfwf, CBF', csfwf, CSF', index_range, pindex);

% Clear both yyaxes first
yyaxis(hfull.pfwaveform, 'left');
cla(hfull.pfwaveform);
yyaxis(hfull.pfwaveform, 'right');
cla(hfull.pfwaveform);

if METADATA.PATCH % Standard deviation patch ON 

    yyaxis(hfull.pfwaveform, 'left');
    [X, Y] = calc_patch(WFPS.CSF.all, cardiacCycle);
    plot_waveform_patched(hfull.pfwaveform, X, Y, CSF, cardiacCycle, METADATA.c1, METADATA.PLW);
    ylabel(ax, 'Flow (CSF) (mL/s)', 'FontSize', 16);

    yyaxis(hfull.pfwaveform, 'right');
    [X, Y] = calc_patch(WFPS.CBF.all, cardiacCycle);
    plot_waveform_patched(hfull.pfwaveform, X, Y, CBF, cardiacCycle, METADATA.c2, METADATA.PLW);
    ylabel(ax, 'Flow (CBF) (mL/s)', 'FontSize', 16);

else % Alternatively, just plot median WF and surrounding points 

    yyaxis(hfull.pfwaveform, 'left'); %#ok<*UNRCH>
    plot_waveform_lines(ax, WFPS.CSF, cardiacCycle, METADATA.c1, METADATA.PLW)
    ylabel(ax, 'Flow (CSF) (mL/s)', 'FontSize', 16);

    yyaxis(hfull.pfwaveform, 'right');
    plot_waveform_lines(ax, WFPS.CBF, cardiacCycle, METADATA.c2, METADATA.PLW)
    ylabel(ax, 'Flow (CBF) (mL/s)', 'FontSize', 16);

end

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

% Labels for centerline visualization? 
[branchLabeled, Ntxt, txt] = calc_labels(branchLabeled, bnum, branchList, Labeltxt, value, average);

% --- Executes when user attempts to close ParameterTool.
function ParameterTool_CloseRequestFcn(hObject, ~, ~)
delete(hObject);

