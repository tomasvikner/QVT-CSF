function varargout = paramMap(varargin)
% PARAMMAP MATLAB code for paramMap.fig
%   Allows for 3D interaction and visualization of 4D flow MRI-derived
%   hemodynamic parameters.
%
%   PARAMMAP, by itself, creates a new PARAMMAP or raises the existing
%   singleton*.
%
%   H = PARAMMAP returns the handle to a new PARAMMAP or the handle to
%   the existing singleton*.
%
%   PARAMMAP('CALLBACK',hObject,eventData,handles,...) calls the local
%   function named CALLBACK in PARAMMAP.M with the given input arguments.
%
%   PARAMMAP('Property','Value',...) creates a new PARAMMAP or raises the
%   existing singleton*.  Starting from the left, property value pairs are
%   applied to the GUI before paramMap_OpeningFcn gets called.  An
%   unrecognized property name or invalid value makes property applicationLeft ICA Cavernous
%   stop.  All inputs are passed to paramMap_OpeningFcn via varargin.
%
%   *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%   instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES
% Edit the above text to modify the response to help paramMap
% Last Modified by GUIDE v2.5 05-Feb-2022 08:20:40

% Developed by Carson Hoffman and Grant Roberts
% University of Wisconsin-Madison 2019
%   Used by: NONE (START FILE)
%   Dependencies: loadpcvipr.m

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
end
% End initialization code - DO NOT EDIT


% --- Executes just before paramMap is made visible.
function paramMap_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to paramMap (see VARARGIN)

disp('paramMap QVT-CSF'); 
disp('Reminder to load case: '); 

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
%set(handles.ParameterTool,'Position',[81 8 190 48]); %HOME
set(handles.TextUpdate,'String','Load in a 4D Flow Dataset');

% --- Executes on button press in LoadData.
function LoadData_Callback(hObject, eventdata, handles)

% disp(['Handles: ' handles])

% hObject    handle to LoadData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Choose default command line output for paramMap
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes paramMap wait for user response (see UIRESUME)
% uiwait(handles.ParameterTool);

% Create global namespace
global branchList Planes hfull p branchLabeled Ntxt nframes res matrix VENC
global directory AveAreaBranch LogPoints fullCData area_val flowPerHeartCycle_val
global PI_val diam_val maxVel_val RI_val flowPulsatile_val timeres segment
global r timeMIPcrossection segmentFull vTimeFrameave velMean_val versionNum
global dcm_obj fig hpatch hscatter Labeltxt cbar hDataTip SavePath
global MAGcrossection bnumMeanFlow bnumStdvFlow StdvFromMean
global VplanesAllx VplanesAlly VplanesAllz imageData caseFilePath
global vesselsAnalyzed allNotes
global VplanesCSF flowCSF StructCS CSFSEG CC CSFROI % add CSF 
global FileNameFlow %  
global venc1 venc2 % 
global wfState % change between CSF/CUBE WFs etc 
global WFPS % waveform point save struct
global T % tangents 
global MAXLAG 
MAXLAG = 250; % max CBF to CSF lag (500 when restricted to >0?)
wfState = 'CnoB'; % set initial CSF WF state 
venc1 = 1; venc2 = 1; % should be set earlier 
FileNameFlow = 'C.h5';

% Initial Variables
hfull = handles;
versionNum = 'v1-2'; %paramMap Version
branchLabeled = 0; %used in cursor updatefunction
Ntxt = []; %used in cursor updatefunction
p = []; %used in cursor updatefunction
directory = uigetdir; %interactive directory selection
vesselsAnalyzed = {};
allNotes = cell(length(get(handles.NamePoint,'String')),1);

% Creates list of all .mat files in selected directory
d = dir([directory filesep '*.mat']);
fn = [{'Load New Case'},{d.name}];
[fileIndx,~] = listdlg('PromptString','Select a file:', ...
    'ListSize',[200 300],'SelectionMode','single','ListString',fn);
disp(['fileIndx: ' num2str(fileIndx)])

%%% Data Loading
if  fileIndx > 1  %if a pre-processed case is selected
    
    set(handles.TextUpdate,'String','Loading Preprocessed Data'); drawnow;
    caseFilePath = [directory filesep fn{fileIndx}];
    
    % separately load larger than usual data struct 
    data_structFilePath = [caseFilePath(1:end-4), '-ds.mat'];
    load(data_structFilePath,'data_struct') %load data_struct (save may req edits from QVT to -v7.3 save)
    load(caseFilePath,'Vel_Time_Res') %load data_struct

    % This will be the name used for the Excel file
    finalFolder = regexp(directory,filesep,'split');
    SummaryName = fn{fileIndx};
    SummaryName = [finalFolder{end} '_' SummaryName(1:end-4)];

    % Makes directory if it does already exist (folder is time-stamped)
    warning off
    mkdir(directory,SummaryName);
    % SavePath = [directory filesep SummaryName];
    SavePath = [directory '/PointsV4/']; % edit 
    disp(SavePath); % TVR

    % Create excel files save summary data
    col_header = ({'Vessel Label', 'Centerline Point', 'Notes',['Max Velocity < ' num2str(VENC) 'cm/s'], ...
        'Mean Flow ml/s','Pulsatility Index','Branch Number'});
    % xlwrite([SavePath filesep 'SummaryParamTool.xls'],col_header,'Summary_Centerline','A1');
    % xlwrite([SavePath filesep 'SummaryParamTool.xls'],get(handles.NamePoint,'String'),'Summary_Centerline','A2');
    
    % New Data Structure
    area_val = data_struct.area_val; %area of vessels
    diam_val = data_struct.diam_val; %diameter of vessels
    branchList = data_struct.branchList; %point locations/labelings
    flowPerHeartCycle_val = data_struct.flowPerHeartCycle_val; %TA flow
    maxVel_val = data_struct.maxVel_val; %TA max velocities
    velMean_val = data_struct.velMean_val; %TA mean velocities
    nframes = data_struct.nframes; %number of temporal cardiac frames
    matrix = data_struct.matrix; %image matrix size (pixels)
    res = data_struct.res; %image resolution (mm)
    timeres = data_struct.timeres; %temporal resolution (ms)
    VENC = data_struct.VENC;
    segment = data_struct.segment; %binary mask (angiogram)
    PI_val = data_struct.PI_val; %pulsatility index
    RI_val = data_struct.RI_val; %resistivity index
    flowPulsatile_val = data_struct.flowPulsatile_val; %TR flow
    r = data_struct.r; %radius of plane (plane size=(2*r)+1)
    timeMIPcrossection = data_struct.timeMIPcrossection; %complex diff.
    MAGcrossection = data_struct.MAGcrossection; %magnitude (in-plane)
    segmentFull = data_struct.segmentFull; %cross-sectional plane masks
    vTimeFrameave = data_struct.vTimeFrameave; %velocity (in-plane)
    Planes = data_struct.Planes; %outer coordinates of plane
    bnumMeanFlow = data_struct.bnumMeanFlow; %mean flow along branches
    bnumStdvFlow = data_struct.bnumStdvFlow; %stdv flow of branches
    StdvFromMean = data_struct.StdvFromMean; %CoV along branches

    VplanesAllx = Vel_Time_Res.VplanesAllx; %TR vel planes (uninterped)
    VplanesAlly = Vel_Time_Res.VplanesAlly;
    VplanesAllz = Vel_Time_Res.VplanesAllz;

    % add for CSF 
    VplanesCSF = []; 
    StructCS = data_struct.StructCS; 
    VplanesCSF.x = Vel_Time_Res.VplanesCSF.x; %TR vel planes (uninterped)
    VplanesCSF.y = Vel_Time_Res.VplanesCSF.y;
    VplanesCSF.z = Vel_Time_Res.VplanesCSF.z;
    fns = fieldnames(data_struct.flowCSF);
    for fi = 1:numel(fns)
        fn = fns{fi};
        flowCSF.(fn).median = data_struct.flowCSF.(fn).median;
        flowCSF.(fn).mean = data_struct.flowCSF.(fn).mean;
    end
    CSFSEG = data_struct.CSFSEG;
    CSFROI = data_struct.CSFROI;
    
    set(handles.TextUpdate,'String','Loading Complete'); drawnow;
    pause(1)
    set(handles.TextUpdate,'String','Please Select Analysis Plane Location'); drawnow;

else % Load in pcvipr data from scratch
    VplanesCSF = []; 
    if exist([directory filesep FileNameFlow],'file') % 
        disp('Loading HDF5')
        [nframes,matrix,res,timeres,VENC,area_val,diam_val,flowPerHeartCycle_val, ...
        maxVel_val,PI_val,RI_val,flowPulsatile_val,velMean_val, ...
        VplanesAllx,VplanesAlly,VplanesAllz,Planes,branchList,segment,r, ...
        timeMIPcrossection,segmentFull,vTimeFrameave,MAGcrossection, imageData, ...
        bnumMeanFlow,bnumStdvFlow,StdvFromMean, ...
        VplanesCSF, flowCSF, StructCS, CSFSEG, T, CSFROI] ...
        = loadHDF5(directory,handles,FileNameFlow);
        % = loadHDF5_py(directory,handles); 
    elseif exist([directory filesep 'CD.dat'],'file')
        [nframes,matrix,res,timeres,VENC,area_val,diam_val,flowPerHeartCycle_val, ...
        maxVel_val,PI_val,RI_val,flowPulsatile_val,velMean_val, ...
        VplanesAllx,VplanesAlly,VplanesAllz,Planes,branchList,segment,r, ...
        timeMIPcrossection,segmentFull,vTimeFrameave,MAGcrossection, imageData, ...
        bnumMeanFlow,bnumStdvFlow,StdvFromMean, T, CSFROI] ...  
        = loadpcvipr(directory,handles); 
    end 
    
    % directory = uigetdir; %select saving dir % load/save in same
    % Save all variables needed to run parametertool. This will be used
    % later to load in data faster instead of having to reload all data.
    % Save data_structure with time/version-stamped filename in 'directory'
    time = datestr(now);
    saveState = [time(1:2) time(4:6) time(10:11) '_' time(13:14) time(16:17)];
    % saveState = [time(1:2) time(4:6) time(10:11) '_' time(13:14) time(16:17) '_' versionNum];
    set(handles.TextUpdate,'String',['Saving Data as Flow4D_' saveState '.mat']); drawnow;

    disp('Assigning values to data_struct...'); 
    
    data_struct = [];
    data_struct.directory = directory;
    data_struct.area_val = area_val;
    data_struct.diam_val = diam_val;
    data_struct.branchList = branchList;
    data_struct.flowPerHeartCycle_val = flowPerHeartCycle_val;
    data_struct.maxVel_val = maxVel_val;
    data_struct.velMean_val = velMean_val;
    data_struct.nframes = nframes;
    data_struct.matrix = matrix;
    data_struct.res = res;
    data_struct.timeres = timeres;
    data_struct.VENC = VENC;
    data_struct.segment = segment;
    data_struct.PI_val = PI_val;
    data_struct.RI_val = RI_val;
    data_struct.flowPulsatile_val = flowPulsatile_val;
    
    data_struct.r = r;
    data_struct.timeMIPcrossection = timeMIPcrossection;
    data_struct.MAGcrossection = MAGcrossection;
    data_struct.segmentFull = segmentFull;
    data_struct.vTimeFrameave = vTimeFrameave;
    data_struct.Planes = Planes;
    data_struct.bnumMeanFlow = bnumMeanFlow;
    data_struct.bnumStdvFlow = bnumStdvFlow;
    data_struct.StdvFromMean = StdvFromMean;
    Vel_Time_Res.VplanesAllx = VplanesAllx; %TR vel planes (uninterped)
    Vel_Time_Res.VplanesAlly = VplanesAlly;
    Vel_Time_Res.VplanesAllz = VplanesAllz;

    % add for CSF
    data_struct.StructCS = StructCS; 
    Vel_Time_Res.VplanesCSF.x = VplanesCSF.x; %TR vel planes (uninterped)
    Vel_Time_Res.VplanesCSF.y = VplanesCSF.y;
    Vel_Time_Res.VplanesCSF.z = VplanesCSF.z;
    fns = fieldnames(flowCSF);
    for fi = 1:numel(fns)
        fn = fns{fi};
        data_struct.flowCSF.(fn).median = flowCSF.(fn).median;
        data_struct.flowCSF.(fn).mean = flowCSF.(fn).mean;
    end
    data_struct.CSFSEG = CSFSEG; 
    data_struct.CSFROI = CSFROI; 
    
    % Saves processed data in same location as pcvipr.mat files
    caseFilePath = fullfile(directory,['Flow4D_' saveState '.mat']);

    % removing data_struct to its own save due to size limits
    % save(caseFilePath,'data_struct','Vel_Time_Res','imageData')
    save(caseFilePath,'Vel_Time_Res','imageData')
    data_structFilePath = [caseFilePath(1:end-4), '-ds.mat'];
    save(data_structFilePath, 'data_struct', '-v7.3') % if too large  
    
    % This will be the name used for the Excel file
    finalFolder = regexp(directory,filesep,'split');
    SummaryName = [finalFolder{end} '_Flow4D_' saveState];
    warning off

    % Where to save data images and excel summary files
    SavePath = [directory filesep SummaryName];   

end

%%% Plotting 3D Interactive Display
set(handles.parameter_choice,'Value',3); %set parameter to flow as default
set(handles.Transparent, 'Value',0);
set(handles.AreaThreshSlide, 'Value',0);

% Initialize visualization
fig = figure(1); cla
%set(fig,'Position',[2325 57 1508 1047]); %WORK
%set(fig,'Position',[1856 37 1416 954]); %HOME

hpatch = patch(isosurface(permute(segment,[2 1 3]),0.5),'FaceAlpha',0); %bw iso angiogram
reducepatch(hpatch,0.7);
set(hpatch,'FaceColor','white','EdgeColor', 'none','PickableParts','none');
set(gcf,'color','black');
axis off tight
% view([1 0 0]);
view([0 0 1]);
axis vis3d
daspect([1 1 1])
set(gca,'zdir','reverse')
camlight headlight;
lighting gouraud
colorbar('off')

% Turn on data cursormode within the figure
dcm_obj = datacursormode(fig); %create dataCursorManager object
datacursormode on;
dcm_obj.DisplayStyle = 'window';
set(handles.CBARmin,'String','min')
set(handles.CBARmax,'String','max')
branchLabeled = 0;

% This will be used in the update function for cursor text
Labeltxt = {'Flow: ',  ' mL/s ';'Average: ',' mL/s '};
cdata = flowPerHeartCycle_val;
hold on
dotSize = 25;
hscatter = scatter3(branchList(:,1),branchList(:,2),branchList(:,3),dotSize,cdata,'filled');
hold off

caxis([min(cdata) max(cdata)]);
cbar = colorbar;
caxis([0 0.8*max(flowPerHeartCycle_val(:))])
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
    'Parent', fig.CurrentAxes); %fill3(pty',ptx',ptz','r') for isosurface
hold off

% Update string (undocumentedmatlab.com/articles/controlling-plot-data-tips)
set(dcm_obj,'UpdateFcn',@myupdatefcn_all); %update dataCursor w/ cust. fcn
hDataTip = dcm_obj.createDatatip(hscatter);

% Convert toolbar to old style and add hot keys
fig.CurrentAxes.Toolbar = [];
addToolbarExplorationButtons(fig)
updateDataCursors(dcm_obj)

% Calculate average area per branch
AveAreaBranch = size(max(branchList(:,4)),1);
for n = 1:max(branchList(:,4))
    Btemp = branchList(:,4)==n;
    AveAreaBranch(n,1) = mean(area_val(Btemp)); %mean area of branch
end

AveAreaBranch = size(max(branchList(:,4)),1);
for n = 1:max(branchList(:,4))
    Btemp = branchList(:,4)==n;
    AveAreaBranch(n,1) = mean(area_val(Btemp)); %mean area of branch
end

LogPoints = true(size(branchList,1),1); %logical array of 1s for areaThresh
fullCData = flowPerHeartCycle_val; %initialize fullCData color as flow

steps = [1./(nframes-1) 10./(nframes-1)]; %set so one 'slide' moves to the next slice exactly
set(handles.VcrossTRslider,'SliderStep',steps);

% new initializations to track manual coregis and segmentations 
nplanes = size(branchList, 1);
CSFSEG.coregTrack = zeros(nplanes, 1);
CSFSEG.mansegTrack = zeros(nplanes, 1);
wfState = 'PC-1'; % start at PCA waveform
disp(['LoadData_Callback complete; WF state: ' wfState])

% --- Outputs from this function are returned to the command line.
function varargout = paramMap_OutputFcn(hObject, eventdata, handles)
varargout{1} = handles.output;


% --- Executes on selection change in parameter_choice.
function parameter_choice_Callback(hObject, eventdata, handles)
global hscatter area_val RI_val PI_val dcm_obj fig velMean_val
global flowPerHeartCycle_val fullCData maxVel_val cbar Labeltxt diam_val
global StdvFromMean
global CC_val % add coupling 

disp('Called parameter_choice_Callback')

% Get parameter option
val = get(handles.parameter_choice, 'Value');
str = get(handles.parameter_choice, 'String');
switch str{val}
    case 'Area'
        % This will be used in the update function for cursor text
        Labeltxt = {'Area: ',  ' cm^2';'Average: ',' cm^2'};
        hscatter.CData = area_val; %update colors on centerline display
        caxis(fig.CurrentAxes,[0 1.5*mean(hscatter.CData)])
        cl = caxis(fig.CurrentAxes);
        set(get(cbar,'xlabel'),'string','Area (cm^2)','fontsize',16,'Color','white');
        set(cbar,'FontSize',16,'color','white');
        fullCData = area_val;
    case 'Ratio of Areas'  
        Labeltxt = {'Area Ratio: ',  ' ';'Average: ',' '};
        disp('Setting Ratio of Areas to CC_val'); 
        % hscatter.CData = diam_val;
        hscatter.CData = CC_val; % 
        caxis(fig.CurrentAxes,[min(hscatter.CData) max(hscatter.CData)]);
        cl = caxis(fig.CurrentAxes);
        % set(get(cbar,'xlabel'),'string','Area Ratio','fontsize',16,'Color','white');
        set(get(cbar,'xlabel'),'string','CSF Coupling','fontsize',16,'Color','white'); 
        set(cbar,'FontSize',16,'color','white');
        % fullCData = diam_val; 
        fullCData = CC_val; 
    case 'Total Flow'
        Labeltxt = {'Flow: ',  ' mL/s';'Average: ',' mL/s'};
        hscatter.CData = flowPerHeartCycle_val;
        caxis(fig.CurrentAxes,[0 0.8*max(flowPerHeartCycle_val(:))])
        cl = caxis(fig.CurrentAxes);
        set(get(cbar,'xlabel'),'string','Flow (mL/s)','fontsize',16,'Color','white');
        set(cbar,'FontSize',16,'color','white');
        fullCData = flowPerHeartCycle_val;
    case 'Maximum Velocity '
        Labeltxt = {'Max Velocity: ',  ' cm/s';'Average: ',' cm/s'};
        hscatter.CData = maxVel_val;
        caxis(fig.CurrentAxes,[min(hscatter.CData) max(hscatter.CData)])
        cl = caxis(fig.CurrentAxes);
        set(get(cbar,'xlabel'),'string','Max Velocity (cm/s)','fontsize',16,'Color','white');
        set(cbar,'FontSize',16,'color','white');
        fullCData = maxVel_val;
    case 'Mean Velocity'
        Labeltxt = {'Mean Velocity: ',  ' cm/s';'Average: ',' cm/s'};
        hscatter.CData = velMean_val;
        caxis(fig.CurrentAxes,[min(hscatter.CData) max(hscatter.CData)])
        cl = caxis(fig.CurrentAxes);
        set(get(cbar,'xlabel'),'string','Mean Velocity (cm/s)','fontsize',16,'Color','white');
        set(cbar,'FontSize',16,'color','white');
        fullCData = velMean_val;
    case 'Flow Consistency'
        Labeltxt = {'Flow Consistency Metric: ',  ' ';'Stdv from Mean: ',' '};
        hscatter.CData = StdvFromMean;
        caxis(fig.CurrentAxes,[0 4])
        cl = caxis(fig.CurrentAxes);
        set(get(cbar,'xlabel'),'string','Flow Consistency Metric','fontsize',16,'Color','white');
        set(cbar,'FontSize',16,'color','white');
        fullCData = StdvFromMean;
    case 'Resistance Index'
        Labeltxt = {'Resistance Index: ',  ' ';'Average: ',' '};
        hscatter.CData = RI_val;
        caxis(fig.CurrentAxes,[-0.5 1])
        cl = caxis(fig.CurrentAxes);
        set(get(cbar,'xlabel'),'string','Resistance Index','fontsize',16,'Color','white');
        set(cbar,'FontSize',16,'color','white');
        fullCData = RI_val;
    case 'CSF Coupling: '
        Labeltxt = {'CSF Coupling: ',  ' ';'Average: ',' '};
        hscatter.CData = CC_val;
        caxis(fig.CurrentAxes,[-1 1])
        cl = caxis(fig.CurrentAxes);
        set(get(cbar,'xlabel'),'string','CSF Coupling','fontsize',16,'Color','white');
        set(cbar,'FontSize',16,'color','white');
        fullCData = CC_val;
    case str{val}
        Labeltxt = {'Pulsatility Index: ',  ' ';'Average: ',' '};
        hscatter.CData = PI_val;
        caxis(fig.CurrentAxes,[0 2])
        cl = caxis(fig.CurrentAxes);
        set(get(cbar,'xlabel'),'string','Pulsatility Index','fontsize',16,'Color','white');
        set(cbar,'FontSize',16,'color','white');
        fullCData = PI_val;
end

set(handles.CBARmin,'String',num2str(cl(1)))
set(handles.CBARmax,'String',num2str(cl(2)))
set(handles.CBARmin,'Value',cl(1))
set(handles.CBARmax,'Value',cl(2))
updateDataCursors(dcm_obj)

% --- Executes during object creation, after setting all properties.
function parameter_choice_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function plot_flowWaveform_CreateFcn(hObject, eventdata, handles) % note this function is never called and does nothing? 
disp('Called plot_flowWaveform_CreateFcn; doing nothing?')

% --- Executes on slider movement.
function Transparent_Callback(hObject, eventdata, handles)
global hpatch Sval

Sval = get(hObject,'Value');
set(hpatch,'FaceAlpha',Sval);

% --- Executes during object creation, after setting all properties.
function Transparent_CreateFcn(hObject, eventdata, handles)
global Sval

if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
set(hObject, 'Min', 0);
set(hObject, 'Max', 1);
set(hObject,'Value',0);
Sval = get(hObject,'Value');


function CBARmin_Callback(hObject, eventdata, handles)
global fig

maxV = str2double(get(handles.CBARmax,'String'));
minV = str2double(get(handles.CBARmin,'String'));
caxis(fig.CurrentAxes,[minV maxV])

% --- Executes during object creation, after setting all properties.
function CBARmin_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function CBARmax_Callback(hObject, eventdata, handles)
global fig

maxV =   str2double(get(handles.CBARmax,'String'));
minV =   str2double(get(handles.CBARmin,'String'));
caxis(fig.CurrentAxes,[minV maxV])

% --- Executes during object creation, after setting all properties.
function CBARmax_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in CBARselection.
function CBARselection_Callback(hObject, eventdata, handles)
global fig

contents = cellstr(get(hObject,'String')); %turn color options to cells
colormap(fig.Children(end),contents{get(hObject,'Value')})

% --- Executes during object creation, after setting all properties.
function CBARselection_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in NamePoint.
function NamePoint_Callback(hObject, eventdata, handles)
global PointLabel dcm_obj vesselsAnalyzed

contents = cellstr(get(handles.NamePoint,'String'));
PointLabel = contents{get(handles.NamePoint,'Value')};
if sum(contains(vesselsAnalyzed,PointLabel))
    set(handles.NamePoint,'ForegroundColor',[0.6 0.6 0.6]);
else
    set(handles.NamePoint,'ForegroundColor',[0 0 0]);
end 
updateDataCursors(dcm_obj)

% --- Executes during object creation, after setting all properties.
function NamePoint_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
global PointLabel

contents = cellstr(get(hObject,'String'));
PointLabel = contents{get(hObject,'Value')};

function toggleCnoB_Callback(hObject, eventdata, handles)
global wfState dcm_obj
if ~strcmp(wfState, 'CnoB')
    disp('WF switch to CnoB');
    wfState = 'CnoB';
end
set(dcm_obj,'UpdateFcn',@myupdatefcn_all); 

function togglePCA_Callback(hObject, eventdata, handles)
global wfState dcm_obj
if ~strcmp(wfState, 'PC-1')
    disp('WF switch to PC-1');
    wfState = 'PC-1';
end
set(dcm_obj,'UpdateFcn',@myupdatefcn_all); 
 
function toggleCUBE_Callback(hObject, eventdata, handles)
global wfState dcm_obj
if ~strcmp(wfState, 'CUBE')
    disp('WF switch to CUBE');
    wfState = 'CUBE';
end
set(dcm_obj,'UpdateFcn',@myupdatefcn_all); 

function updateWaveforms(ROITYPE)
global VplanesCSF CSFSEG dcm_obj r

nframes = 20; 
[INDEX, ~] = get_pindex(dcm_obj);

% Extract and reshape the ROI mask
roiCSF = CSFSEG.(ROITYPE)(INDEX, :);
imdim = sqrt(numel(roiCSF));
roiCSF = reshape(roiCSF, imdim, imdim);

% Get velocity plane field names, excluding tracking fields
fns = setdiff(fieldnames(CSFSEG), {'coregTrack', 'segTrack', 'mthrTrack', 'sthrTrack', 'mansegTrack'}, 'stable');

% Loop over time frames
for n = 1:nframes
    % Extract and resize velocity planes
    v = struct();
    for axis = ["x", "y", "z"]
        fld = char(axis);
        slice = squeeze(VplanesCSF.(fld)(INDEX, :, n));             % Extract 1D vector
        reshaped = reshape(slice, (2 * r) + 1, []);                  % Reshape to 2D
        resized = imresize(reshaped, [imdim imdim]);                % Resize to image dimension
        v.(fld) = resized;
    end

    % Combine into total velocity
    vtotal = 0.1 * (v.x + v.y + v.z);

    % Apply CSF mask
    maskedV = roiCSF(:) .* vtotal(:);

    % Compute flow stats for each ROI
    for i = 1:numel(fns)
        fn = fns{i};
        valid = maskedV(maskedV ~= 0);  % Use only non-zero (i.e., masked-in) values
        flowCSF.(fn).median(INDEX, n) = median(valid, 'omitnan');
        flowCSF.(fn).mean(INDEX, n)   = mean(valid, 'omitnan');
    end
end

function manualCoreg_Callback(~, ~, ~)
global dcm_obj CSFSEG segmentFull branchList

disp('Called manual 2D point-point Co-Reg')
pindex = getBranchIndices(dcm_obj, branchList);
imdim = sqrt(numel(segmentFull(1, :)));

moving = reshape(CSFSEG.cube(pindex, :), imdim, imdim);
fixed  = reshape(CSFSEG.mcsf(pindex, :), imdim, imdim); % Using mag only

figure(11); imshow(moving, []); [xm, ym] = ginput(1);
figure(12); imshow(fixed,  []); [xf, yf] = ginput(1);

dx = xf - xm; dy = yf - ym;
translated = imtranslate(moving, [dx, dy], 'OutputView', 'same');
figure(13); imshow(translated, []);

CSFSEG.cube(pindex, :)  = translated(:);
CSFSEG.bcube(pindex, :) = translated(:) > graythresh(translated(:));
CSFSEG.coregTrack(pindex) = 1;

set(dcm_obj, 'UpdateFcn', @myupdatefcn_all);
close([11 12 13]);
updateWaveforms('dcube');

function pindex = getBranchIndices(dcm_obj, branchList)
info = getCursorInfo(dcm_obj);
ptList = reshape([info.Position], 3, []).';
pindex = zeros(size(ptList, 1), 1);

for n = 1:size(ptList,1)
    idx = all(bsxfun(@eq, branchList(:,1:3), ptList(n,:)), 2);
    pindex(n) = find(idx, 1);
end

function manualCoreg_Callback_3D(~, ~, ~)
global dcm_obj CSFSEG segmentFull branchList

disp('Called manual 3D point-point Co-Reg')
pindex = getBranchIndices(dcm_obj, branchList);
imdim = sqrt(numel(segmentFull(1,:)));

moving = reshape(CSFSEG.cube(pindex,:), imdim, imdim);
fixed  = reshape(CSFSEG.mcsf(pindex,:), imdim, imdim);

% Get multiple matching points
figure(11); imshow(moving, []); [xm, ym] = ginput(5);
figure(12); imshow(fixed,  []); [xf, yf] = ginput(5);

% Estimate and apply rigid 2D transformation
tform = fitgeotrans([xm ym], [xf yf], 'nonreflectivesimilarity');
translated = imwarp(moving, tform, 'OutputView', imref2d(size(moving)));

% Display result
figure(13); imshow(translated, []);

% Update segmentation
CSFSEG.cube(pindex,:)  = translated(:);
CSFSEG.bcube(pindex,:) = translated(:) > graythresh(translated(:));
CSFSEG.coregTrack(pindex) = 1;

set(dcm_obj, 'UpdateFcn', @myupdatefcn_all);
close([11 12 13]);
updateWaveforms('dcube');

function manualSeg_Callback(~, ~, ~)
global dcm_obj CSFSEG segmentFull branchList

disp('Called manual 2D segmentation')
pindex = getBranchIndices(dcm_obj, branchList);
imdim = sqrt(numel(segmentFull(1,:)));

% moving = reshape(CSFSEG.cube(pindex,:), imdim, imdim);
moving = reshape(CSFSEG.scsf(pindex,:), imdim, imdim); % Manual seg on cnob, not cube 

figure(11); imshow(moving, []);
disp('Draw region with freehand tool. Double-click to finish.');
caxis([min(moving(:)), 0.5 * max(moving(:))]);

h = imfreehand;
mask = createMask(h);
close(11);

% Update segmentation
% CSFSEG.bcube(pindex,:) = mask(:);
CSFSEG.cnob(pindex,:) = mask(:); % Manual seg on cnob, not cube
CSFSEG.segTrack(pindex) = 1;

set(dcm_obj, 'UpdateFcn', @myupdatefcn_all);
updateWaveforms('cnob');

% --- Executes on button press in SavePoint.
function SavePoint_Callback(hObject, eventdata, handles)
global PointLabel nframes VENC timeres branchList timeMIPcrossection area_val
global flowPerHeartCycle_val PI_val diam_val maxVel_val RI_val flowPulsatile_val
global vTimeFrameave velMean_val dcm_obj fig segmentFull SavePath MAGcrossection
global vesselsAnalyzed allNotes 
global venc1 venc2 flowCSF directory % 
global WFPS VcrossTR CcrossTR Maskcross CC % waveform point save struct
% global CC_val % CBF-CSF coupling coefficient (not in use) 

% CSFROI values need only be called when using the save callback
global CSFROI

vesselsAnalyzed{end+1} = PointLabel;

set(handles.TextUpdate,'String','Saving Data.');drawnow;

info_struct = getCursorInfo(dcm_obj);
ptList = [info_struct.Position];
ptList = reshape(ptList,[3,numel(ptList)/3])';
pindex = zeros(size(ptList,1),1);

for n = 1:size(ptList,1)
    xIdx = find(branchList(:,1) == ptList(n,1));
    yIdx = find(branchList(xIdx,2) == ptList(n,2));
    zIdx = find(branchList(xIdx(yIdx),3) == ptList(n,3));
    pindex(n) = xIdx(yIdx(zIdx));
end

% Gives associated branch number if full branch point is wanted
bnum = branchList(pindex,4);
Logical_branch = branchList(:,4) ~= bnum;

% OUTPUT +/- points use this
index_range = pindex-2:pindex+2;

%removes outliers and points from other branches
index_range(index_range<1) = [];
index_range(index_range>size(branchList,1)) = [];
index_range(Logical_branch(index_range)) = [];

% disp('flowPulsatile (paramMap.m)')
% Time-resolved flow
flowPulsatile = flowPulsatile_val(index_range,:);
flowPulsatile = [flowPulsatile;mean(flowPulsatile,1);std(flowPulsatile,1)];

% savepoint_callback add CSF flow
NOXLS = 1; 
if NOXLS 

    savePoint = [];

    savePoint.TRCBF = VcrossTR;
    savePoint.TRCSF = CcrossTR;
    % savePoint.BMASK = Maskcross;
    savePoint.CMASK = CC.madj;
    
    savePoint.bnum = bnum;
    savePoint.Logical_branch = Logical_branch;
    savePoint.index_range = index_range; 
    savePoint.ptList = ptList; 

    savePoint.CBF = flowPulsatile * venc2; 
    savePoint.CSF.bics = flowCSF.bics.mean(index_range, :) * venc1; 
    savePoint.CSF.bcsf = flowCSF.bcsf.mean(index_range, :) * venc1; 
    savePoint.CSF.cnob = flowCSF.cnob.mean(index_range, :) * venc1; 
    savePoint.CSF.full = flowCSF.full.mean(index_range, :) * venc1; 
    savePoint.CSF.pc1 = flowCSF.PC1.mean(index_range, :) * venc1; 
    DOMADJ = false; 
    if DOMADJ
        if sum(flowCSF.madj.mean(index_range, :)) > 0 %#ok<UNRCH> 
            savePoint.CSF.madj = flowCSF.madj.mean(index_range, :) * venc1; 
        end
        savePoint.CSF.madj = flowCSF.madj.mean(index_range, :) * venc1; 
    end

    savePoint.CSFROI.bics = CSFROI.bics(index_range, :, :, :); 
    savePoint.CSFROI.bcsf = CSFROI.bcsf(index_range, :, :, :); 
    savePoint.CSFROI.cnob = CSFROI.cnob(index_range, :, :, :); 
    savePoint.CSFROI.flow = CSFROI.flow(index_range, :, :, :); 

    % Name of point label e.g. Left MCA as in scroll menu in GUI 
    savePoint.pointName = PointLabel; 

    % note - Waveform point save (struct) - THIS contains most important information to export!
    savePoint.WFPS = WFPS;  

    % save([SavePath savePoint.pointName '.mat'], 'savePoint');
    SavePath = [directory filesep 'PointsV4/']; 
    if ~exist(SavePath, 'dir')
        mkdir(SavePath)
    end
    save([SavePath savePoint.pointName '.mat'], 'savePoint');
    disp(['Saved point @: ' [SavePath savePoint.pointName '.mat'] ]);

else % temp remove old save to xls file 
    disp('Save to XLS temp cleaned out in CSF V4'); 
end

% --- NoteBox_Callback
function NoteBox_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function NoteBox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in SubmitNote.
function SubmitNote_Callback(hObject, eventdata, handles)
global SavePath allNotes

set(handles.TextUpdate,'String','Saving Note for Current Vessel');drawnow;
SaveLoc =  sprintf('C%i',get(handles.NamePoint,'Value')+1);
Notes = {get(handles.NoteBox,'String')}; %get any notes from notebox
allNotes(get(handles.NamePoint,'Value')+1) = Notes;
xlwrite([SavePath filesep 'SummaryParamTool.xls'],Notes,'Summary_Centerline',SaveLoc);
set(handles.NoteBox,'String',' ');
set(handles.TextUpdate,'String','Done Saving Note');drawnow;

% --- Executes on button press in AxialView.
function AxialView_Callback(hObject, eventdata, handles)
global fig
view(fig.CurrentAxes,[180,90])


% --- Executes on button press in SagittalView.
function SagittalView_Callback(hObject, eventdata, handles)
global fig
view(fig.CurrentAxes,[180,0])


% --- Executes on button press in CoronalView.
function CoronalView_Callback(hObject, eventdata, handles)
global fig
view(fig.CurrentAxes,[90,0])


% --- Executes on slider movement.
function AreaThreshSlide_Callback(hObject, eventdata, handles)
global LogPoints area_val branchList hscatter AveAreaBranch PI_val RI_val
global velMean_val diam_val maxVel_val flowPerHeartCycle_val StdvFromMean
global CC_val %  

LogPoints = find(AveAreaBranch>max(AveAreaBranch)*get(hObject,'Value')*.15);
LogPoints = ismember(branchList(:,4),LogPoints);

if get(handles.InvertArea,'Value') == 0
    hscatter.XData = branchList(LogPoints,1);
    hscatter.YData = branchList(LogPoints,2);
    hscatter.ZData = branchList(LogPoints,3);
    
    val = get(handles.parameter_choice, 'Value');
    str = get(handles.parameter_choice, 'String');
    switch str{val}
        case 'Area'
            hscatter.CData = area_val(LogPoints);
        case 'Ratio of Areas'
            hscatter.CData = CC_val(LogPoints);
            hscatter.CData = CC_val(LogPoints); 
        case 'Total Flow'
            hscatter.CData = flowPerHeartCycle_val(LogPoints);
        case 'Maximum Velocity '
            hscatter.CData = maxVel_val(LogPoints);
        case 'Mean Velocity'
            hscatter.CData = velMean_val(LogPoints);
        case 'Flow Consistency'
            hscatter.CData = StdvFromMean(LogPoints);
        case 'Resistance Index'
            hscatter.CData = RI_val(LogPoints);
        case 'CSF Coupling'
            hscatter.CData = CC_val(LogPoints);
        case str{val}
            hscatter.CData = PI_val(LogPoints);
    end
else
    hscatter.XData = branchList(~LogPoints,1);
    hscatter.YData = branchList(~LogPoints,2);
    hscatter.ZData = branchList(~LogPoints,3);
    
    val = get(handles.parameter_choice, 'Value');
    str = get(handles.parameter_choice, 'String');
    switch str{val}
        case 'Area'
            hscatter.CData = area_val(~LogPoints);
        case 'Ratio of Areas'
            % hscatter.CData = diam_val(~LogPoints);
            hscatter.CData = CC_val(~LogPoints); 
        case 'Total Flow'
            hscatter.CData = flowPerHeartCycle_val(~LogPoints);
        case 'Maximum Velocity '
            hscatter.CData = maxVel_val(~LogPoints);
        case 'Mean Velocity'
            hscatter.CData = velMean_val(~LogPoints);
        case 'Flow Consistency'
            hscatter.CData = StdvFromMean(LogPoints);
        case 'Resistance Index'
            hscatter.CData = RI_val(~LogPoints);
        case 'CSF Coupling'
            hscatter.CData = CC_val(~LogPoints);
        case str{val}
            hscatter.CData = PI_val(~LogPoints);
    end
end


% --- Executes during object creation, after setting all properties.
function AreaThreshSlide_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% --- help function for updating waveforms from sliders/mansegs
function [pindex, index_range] = get_pindex(dcm_obj)
global branchList dcm_obj 
info_struct = getCursorInfo(dcm_obj);
ptList = [info_struct.Position];
ptList = reshape(ptList,[3,numel(ptList)/3])';
pindex = zeros(size(ptList,1),1);
for n = 1:size(ptList,1)
    xIdx = find(branchList(:,1) == ptList(n,1));
    yIdx = find(branchList(xIdx,2) == ptList(n,2));
    zIdx = find(branchList(xIdx(yIdx),3) == ptList(n,3));
    pindex(n) = xIdx(yIdx(zIdx));
end
bnum = branchList(pindex,4);
Logical_branch = branchList(:,4) ~= bnum;
index_range = pindex; 
index_range(index_range<1) = [];
index_range(index_range>size(branchList,1)) = [];
index_range(Logical_branch(index_range)) = [];


% --- Executes on slider movement.
function clWidthSlider_Callback(hObject, eventdata, handles)
global hscatter 
set(hscatter,'SizeData',get(hObject,'Value'));


% --- Executes during object creation, after setting all properties.
function clWidthSlider_CreateFcn(hObject, eventdata, handles)
% disp('Called clWidthSlider_CreateFcn')
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in VisualTool.
function VisualTool_Callback(hObject, eventdata, handles)
global Planes branchList segment caseFilePath res 
set(handles.TextUpdate,'String','Opening Visual Tool'); drawnow;
fourDvis(Planes,branchList,segment,caseFilePath,res);
uiwait;
set(handles.TextUpdate,'String','Visual Tool Closed'); drawnow;


% --- Executes on button press in InvertArea.
function InvertArea_Callback(hObject, eventdata, handles)
% hObject    handle to InvertArea (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of InvertArea
global LogPoints area_val branchList hscatter PI_val RI_val velMean_val
global diam_val maxVel_val flowPerHeartCycle_val StdvFromMean
global CC_val % 

% Capable of inverting areaThresh (keep vessels OUTSIDE/INSIDE areaThresh)
OnOff = get(hObject,'Value'); %on off switch
if OnOff == 0 %if turned off (default),
    hscatter.XData = branchList(LogPoints,1); %plot angio w/in areaThresh
    hscatter.YData = branchList(LogPoints,2);
    hscatter.ZData = branchList(LogPoints,3);
    
    val = get(handles.parameter_choice, 'Value');
    str = get(handles.parameter_choice, 'String');
    switch str{val} %plot centerlines w/in areaThresh
        case 'Area'
            hscatter.CData = area_val(LogPoints);
        case 'Ratio of Areas'
            % hscatter.CData = diam_val(LogPoints); 
            hscatter.CData = CC_val(LogPoints); % 
        case 'Total Flow'
            hscatter.CData = flowPerHeartCycle_val(LogPoints);
        case 'Maximum Velocity '
            hscatter.CData = maxVel_val(LogPoints);
        case 'Mean Velocity'
            hscatter.CData = velMean_val(LogPoints);
        case 'Flow Consistency'
            hscatter.CData = StdvFromMean(LogPoints);
        case 'Resistance Index'
            hscatter.CData = RI_val(LogPoints);
        case 'CSF Coupling'
            hscatter.CData = CC_val(LogPoints);
        case str{val}
            hscatter.CData = PI_val(LogPoints);
    end
else %if invert is turned on, PLOT DATA POINTS OUTSIDE AREA THRESHOLD
    hscatter.XData = branchList(~LogPoints,1);
    hscatter.YData = branchList(~LogPoints,2);
    hscatter.ZData = branchList(~LogPoints,3);
    
    val = get(handles.parameter_choice, 'Value');
    str = get(handles.parameter_choice, 'String');
    switch str{val}
        case 'Area'
            hscatter.CData = area_val(~LogPoints);
        case 'Ratio of Areas'
            hscatter.CData = diam_val(~LogPoints);
        case 'Total Flow'
            hscatter.CData = flowPerHeartCycle_val(~LogPoints);
        case 'Maximum Velocity '
            hscatter.CData = maxVel_val(~LogPoints);
        case 'Mean Velocity'
            hscatter.CData = velMean_val(~LogPoints);
        case 'Flow Consistency'
            hscatter.CData = StdvFromMean(LogPoints);
        case 'Resistance Index'
            hscatter.CData = RI_val(~LogPoints);
        case 'CSF Coupling'
            hscatter.CData = CC_val(~LogPoints);
        case str{val}
            hscatter.CData = PI_val(~LogPoints);
    end
end

% --- Executes on slider movement.
function VcrossTRslider_Callback(hObject, eventdata, handles)
updateVcrossTR(handles)

% --- Executes during object creation, after setting all properties.
function VcrossTRslider_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function updateVcrossTR(handles)
global dcm_obj hfull segmentFull VplanesAllx VplanesAlly VplanesAllz 
global VplanesCSF StructCS CSFSEG CC % add some TR CSF planes 
global nframes branchList wfState

info_struct = getCursorInfo(dcm_obj);
if ~isempty(info_struct)
    ptList = [info_struct.Position];
    ptList = reshape(ptList,[3,numel(ptList)/3])';
    pindex = zeros(size(ptList,1),1);
    % Find cursor point in branchList
    for n = 1:size(ptList,1)
        xIdx = find(branchList(:,1) == ptList(n,1));
        yIdx = find(branchList(xIdx,2) == ptList(n,2));
        zIdx = find(branchList(xIdx(yIdx),3) == ptList(n,3));
        pindex(n) = xIdx(yIdx(zIdx));
    end
    imdim = sqrt(size(segmentFull,2)); %side length of cross-section

    Maskcross = segmentFull(pindex,:);
    Maskcross = reshape(Maskcross,imdim,imdim);

    % get slice number from slider
    sliceNum = 1+round( get(hfull.VcrossTRslider,'Value').*(nframes-1) ); 

    % Orig CS data 
    v1 = squeeze(VplanesAllx(pindex,:,:));
    v2 = squeeze(VplanesAlly(pindex,:,:));
    v3 = squeeze(VplanesAllz(pindex,:,:));
    VcrossTR = 0.1*(v1 + v2 + v3); 

    % CSF CS data
    c1 = squeeze(VplanesCSF.x(pindex,:,:));
    c2 = squeeze(VplanesCSF.y(pindex,:,:));
    c3 = squeeze(VplanesCSF.z(pindex,:,:));
    CcrossTR = 0.1*(c1 + c2 + c3); 

    cbfscale = 1.0; 
    csfscale = 1.0; 
    % VcrossTRmean = mean(VcrossTR, 2); % Time-average
    % VcrossTR = VcrossTR - VcrossTRmean; % remove mean flow to visualize CBF pulsation
    normDim = sqrt(size(VcrossTR,1));
    VcrossTR = reshape(VcrossTR,normDim,normDim,nframes);
    VcrossTR = imresize(VcrossTR,[imdim imdim],'nearest');
    minn = min(Maskcross.*VcrossTR,[],'all')*cbfscale;
    maxx = max(Maskcross.*VcrossTR,[],'all')*cbfscale;
    imshow(VcrossTR(:,:,sliceNum),[minn maxx],'InitialMagnification','fit','Parent',hfull.TRcross)
    title(hfull.TRcross, 'CBF velocity (time-res. )', 'FontSize', 13);
    visboundaries(hfull.TRcross,Maskcross,'LineWidth',1)

    % CSF CS data
    normDim = sqrt(size(CcrossTR,1));
    CcrossTR = reshape(CcrossTR,normDim,normDim,nframes);
    CcrossTR = imresize(CcrossTR,[imdim imdim],'nearest');
    minn = min(CcrossTR(:)); maxx = max(CcrossTR(:));
    minn = minn * csfscale; maxx = maxx * csfscale;
    imshow(CcrossTR(:,:,sliceNum),[minn maxx],'InitialMagnification','fit','Parent',hfull.CSF4)
    title(hfull.CSF4, 'CSF velocity (time-res. )', 'FontSize', 13);
    INCROI = double(StructCS.full(pindex, :)); INCROI = reshape(INCROI,imdim,imdim);
    if strcmp(wfState, 'CUBE')
        visboundaries(hfull.CSF4,CC.bcube,'LineWidth',1) % auto/manual combo to vis/seg
    else
        visboundaries(hfull.CSF4,CC.madj,'LineWidth',1) % auto/manual combo to vis/seg
    end
end 

function txt = myupdatefcn_all(empt,event_obj) %

% closing manual ROI windows when changing pindex 
% disp('myupdatefcn_all')
try 
    close(2:3)
catch
end

% Customizes text of data tips
global Labeltxt branchLabeled PointLabel branchList fullCData
global flowPulsatile_val Planes p dcm_obj Ntxt hfull timeMIPcrossection flowCSF
global segmentFull MAGcrossection vTimeFrameave fig timeres nframes
global VplanesAllx VplanesAlly VplanesAllz
global VplanesCSF StructCS CSFSEG CC % 
global wfState % toggle CSF WF state 
global venc1 venc2 % note these are 1/1 anyway (edit if scaling in QVT)
global WFPS % waveform point save struct
global MAXLAG

info_struct = getCursorInfo(dcm_obj);
ptList = [info_struct.Position];
ptList = reshape(ptList,[3,numel(ptList)/3])';
pindex = zeros(size(ptList,1),1);

% resetMagSlider(hObject)

% Find cursor point in branchList
for n = 1:size(ptList,1)
    xIdx = find(branchList(:,1) == ptList(n,1));
    yIdx = find(branchList(xIdx,2) == ptList(n,2));
    zIdx = find(branchList(xIdx(yIdx),3) == ptList(n,3));
    pindex(n) = xIdx(yIdx(zIdx));
end

% Get associated branch number of full branch
bnum = branchList(pindex,4);
Logical_branch = branchList(:,4) ~= bnum;
index_range = pindex-2:pindex+2; % OUTPUT +/- points
index_range(index_range<1) = []; %removes outliers and other branch points
index_range(index_range>size(branchList,1)) = [];
index_range(Logical_branch(index_range)) = [];

% Update Cross-sectional planes for points
set(p,'XData',Planes(pindex,:,1)','YData',Planes(pindex,:,2)','ZData',Planes(pindex,:,3)')
imdim = sqrt(size(segmentFull,2)); %side length of cross-section

% difining new local CS for vis 
CC = [];
fns = fieldnames(CSFSEG);
for i = 1:numel(fns) 
    fn = fns{i};
    if ~contains(fn, 'Track') % Don't update segmentation and coreg track
        CC.(fn) = reshape(double(CSFSEG.(fn)(pindex,:)),imdim,imdim);
    end
end
if sum(CC.madj(:)==0) 
    CC.madj = CC.auto;
end

% CBF ROI 
Maskcross = segmentFull(pindex,:);
Maskcross = reshape(Maskcross,imdim,imdim); 

% Magnitude TA
MAGcross = MAGcrossection(pindex,:);
MAGcross = reshape(MAGcross,imdim,imdim);
imshow(MAGcross,[],'InitialMagnification','fit','Parent',hfull.MAGcross) % CBF #1 - show CBF mag as normal
title(hfull.MAGcross, 'CBF magnitude', 'FontSize', 13);
% visboundaries(hfull.MAGcross,Maskcross,'LineWidth',1)

% Complex diffference TA
CDcross = timeMIPcrossection(pindex,:);
CDcross = reshape(CDcross,imdim,imdim);
imshow(CC.rage,[],'InitialMagnification','fit','Parent', hfull.CDcross) % CBF #2 - show MPRAGE in CD location 
title(hfull.CDcross, 'T1-w MP-RAGE', 'FontSize', 13);
% visboundaries(hfull.CDcross,Maskcross,'LineWidth',1)

% Velocity TA - through plane
imshow(CDcross,[],'InitialMagnification','fit','Parent', hfull.VELcross) % CBF #3 - show CD in VELcross location 
title(hfull.VELcross, 'Complex Difference', 'FontSize', 13);
visboundaries(hfull.VELcross,Maskcross,'LineWidth',1)

% Velocity TR - through plane
sliceNum = 1+round( get(hfull.VcrossTRslider,'Value').*(nframes-1) ); 
v1 = squeeze(VplanesAllx(pindex,:,:));
v2 = squeeze(VplanesAlly(pindex,:,:));
v3 = squeeze(VplanesAllz(pindex,:,:));
VcrossTR = 0.1*(v1 + v2 + v3); % ORIG 
normDim = sqrt(size(VcrossTR,1));
VcrossTR = reshape(VcrossTR,normDim,normDim,nframes);
VcrossTR = imresize(VcrossTR,[imdim imdim],'nearest');
minn = min(Maskcross.*VcrossTR,[],'all')*1.1;
maxx = max(Maskcross.*VcrossTR,[],'all')*1.1;
imshow(VcrossTR(:,:,sliceNum),[minn maxx],'InitialMagnification','fit','Parent',hfull.TRcross)
title(hfull.TRcross, 'CBF velocity (time-res. )', 'FontSize', 13);
visboundaries(hfull.TRcross,Maskcross,'LineWidth',1) % CBF #4 - show time-resolved CBF 

% CSF vis
imshow(CC.mcsf,[],'InitialMagnification','fit','Parent',hfull.CSF1) % CSF #1 - show CSF mag and blood in CSF (BICS) 
title(hfull.CSF1, 'CSF magnitude', 'FontSize', 13)
visboundaries(hfull.CSF1,CC.bics,'LineWidth',1)
imshow(CC.cube,[],'InitialMagnification','fit','Parent',hfull.CSF2) % CSF #2 - show CUBE and BCSF 
title(hfull.CSF2, 'CUBE Anti-FLAIR', 'FontSize', 13)
% title(hfull.CSF2, 'CUBE Anti-FLAIR', 'FontSize', 13) TODO: what should be
% the final name on the T2 CUBE 
visboundaries(hfull.CSF2,CC.bcube,'LineWidth',1) % is > 0 may be sufficient if interp of non-binary, or fix interp to nn 
imshow(CC.scsf,[],'InitialMagnification','fit','Parent',hfull.CSF3) % CSF #3 - show CSF SD and CNOB 
title(hfull.CSF3, 'CSF velocity STD', 'FontSize', 13)
visboundaries(hfull.CSF3,CC.bcsf,'LineWidth',1)

% CSF Velocity TR - through-plane 
sliceNum = 1+round( get(hfull.VcrossTRslider,'Value').*(nframes-1) ); 
c1 = squeeze(VplanesCSF.x(pindex,:,:));
c2 = squeeze(VplanesCSF.y(pindex,:,:));
c3 = squeeze(VplanesCSF.z(pindex,:,:));
CcrossTR = 0.1*(c1 + c2 + c3); % ORIG tdim = 3; 
normDim = sqrt(size(CcrossTR,1));
CcrossTR = reshape(CcrossTR,normDim,normDim,nframes);
CcrossTR = imresize(CcrossTR,[imdim imdim],'nearest');
minn = min(CcrossTR(:)); maxx = max(CcrossTR(:));
imshow(CcrossTR(:,:,sliceNum),[minn maxx],'InitialMagnification','fit','Parent',hfull.CSF4)
title(hfull.CSF4, 'CSF velocity (time-res. )', 'FontSize', 13)
visboundaries(hfull.CSF4,CC.madj,'LineWidth',1) % CSF #4

% Get value of parameter at point and mean within 5pt window
value = fullCData(pindex);
average = fullCData(index_range);
pside = index_range;
pside(pside==pindex) = [];

% define waveforms 
if size(pside, 2) > size(pside, 1)
    pside = pside';
end
pinds = [pindex; pside];
csfinds = pinds; 
cbfinds = pinds;   

c = linspecer(2); 
c1 = c(1, :); 
c2 = c(2, :);
PLW = 2.4; 

cbfwf = flowPulsatile_val(cbfinds, :);
if strcmp(wfState, 'PC-1')
    csfwf = flowCSF.PC1.mean(csfinds, :); 
elseif strcmp(wfState, 'CnoB')
    csfwf = flowCSF.cnob.mean(csfinds, :); 
elseif strcmp(wfState, 'CUBE')
    csfwf = flowCSF.bcube.mean(csfinds, :); 
end

zeroRows = all(csfwf == 0, 2); % Do not average over those that have no CSF ROI
nanRows = all(isnan(csfwf), 2);
emptyRows = find(zeroRows | nanRows);
csfwf(emptyRows, :) = []; 
cbfwf(emptyRows, :) = []; 

ncsf = size(csfwf, 1);
if ncsf == 0
    CSF = zeros(1, 20);
elseif ncsf == 1
    CSF = csfwf;
else
    % CSF = mean(csfwf);
    CSF = median(csfwf);
end
% CBF = mean(cbfwf);
CBF = median(cbfwf);

nframes = 20; % Original frames
iframes = 1000; % Desired output resolution
xo = linspace(1, nframes, nframes);
xq = linspace(1, nframes, iframes);
CBF = interp1(xo, CBF', xq, 'pchip'); % mean point 
CSF = interp1(xo, CSF', xq, 'pchip');
csfwf = interp1(xo, csfwf', xq, 'pchip'); % all points 
cbfwf = interp1(xo, cbfwf', xq, 'pchip');

[rmax, mlag] = waveformCoupling(CBF, CSF, MAXLAG); % max 300 ms delay?

amp = [];
amp.CSF = max(CSF) - min(CSF);
amp.CBF = max(CBF) - min(CBF);
cdv = [];
ccsf = cumtrapz(CSF - mean(CSF)) * 1e-3; % with current interp? 
ccbf = cumtrapz(CBF - mean(CBF)) * 1e-3;
cdv.CSF = max(ccsf) - min(ccsf);
cdv.CBF = max(ccbf) - min(ccbf);

WFPS = [];
WFPS.amp = amp;
WFPS.cdv = cdv;
WFPS.rmax = rmax;
WFPS.mlag = mlag;
WFPS.CBF.all = cbfwf;
WFPS.CBF.avg = CBF';
WFPS.CSF.all = csfwf;
WFPS.CSF.avg = CSF';
WFPS.irange = index_range;
WFPS.pindex = pindex; 

cardiacCycle = 1:numel(CBF);

if isempty(CSF)
    disp('CSF waveform empty!')
    disp(size(flowCSF.mean))
    disp(size(flowPulsatile_CBF))
end

% Clear both yyaxes first
yyaxis(hfull.pfwaveform, 'left');
cla(hfull.pfwaveform);
yyaxis(hfull.pfwaveform, 'right');
cla(hfull.pfwaveform);

% Now start plotting cleanly
yyaxis(hfull.pfwaveform, 'left');
plot(hfull.pfwaveform, cardiacCycle, WFPS.CSF.avg, 'LineWidth', PLW*2, 'Color', c1);
hold(hfull.pfwaveform, 'on');
for k = 1:size(WFPS.CSF.all, 2)
    plot(hfull.pfwaveform, cardiacCycle, WFPS.CSF.all(:, k), ...
        'LineWidth', PLW/2, 'Color', c1, ...
        'LineStyle', '--', ...
        'Marker', 'none');
end
ylabel(hfull.pfwaveform, 'Flow (CSF) (mL/s)', 'FontSize', 16);
hfull.pfwaveform.YColor = c1;

yyaxis(hfull.pfwaveform, 'right');
plot(hfull.pfwaveform, cardiacCycle, WFPS.CBF.avg, 'LineWidth', PLW*2, 'Color', c2);
hold(hfull.pfwaveform, 'on');
for k = 1:size(WFPS.CBF.all, 2)
    plot(hfull.pfwaveform, cardiacCycle, WFPS.CBF.all(:, k), ...
        'LineWidth', PLW/2, 'Color', c2, ...
        'LineStyle', '-', ...
        'Marker', 'none');
end
ylabel(hfull.pfwaveform, 'Flow (CBF) (mL/s)', 'FontSize', 16);
hfull.pfwaveform.YColor = c2;

% legend(hfull.pfwaveform, '', ...
%     ['CSF amp.: ' num2str(amp.CSF, 3) ' mL/s'], '', ...
%     ['CBF amp.: ' num2str(amp.CBF, 3) ' mL/s'], ...
%     'Box', 'off', 'FontSize', 16, 'FontWeight', 'bold', ...
%     'Location', 'North');

title(hfull.pfwaveform, ttext, 'FontSize', 16);
grid(hfull.pfwaveform, 'on');
hfull.pfwaveform.XAxisLocation = 'bottom';

yyaxis(hfull.pfwaveform, 'left');
axis(hfull.pfwaveform, 'tight');

% Put the number labels on the CenterlinePlot if new branch
if branchLabeled ~= bnum
    delete(Ntxt)
    branchLabeled = bnum;
    index_branch = branchList(:,4) == branchLabeled;
    branchActual = branchList(index_branch,1:3);
    textint = 0:5:length(branchActual)-1;
    numString_val = num2str(textint);
    numString_val = strsplit(numString_val);
    Ntxt = text(branchActual(textint+1,1),branchActual(textint+1,2), ...
        branchActual(textint+1,3),numString_val,'Color','w','FontSize',10,...
        'HitTest','off','PickableParts','none','Parent',fig.Children(2));
end

% Get branch indices and current label point
branchActual = branchList(branchList(:,4) == branchLabeled,5);
CurrentNum = find(branchList(pindex,5)==branchActual)-1;

% Update cursor text
txt = {['Point Label: ' , PointLabel , sprintf('\n'), ...
    Labeltxt{1,1}, sprintf('%0.3f',value),Labeltxt{1,2}, sprintf('\n'), ...
    Labeltxt{2,1},sprintf('%0.3f',mean(average)),Labeltxt{2,2},sprintf('\n'), ...
    'Current Branch #: ',sprintf('%i',CurrentNum),sprintf('\n') ...
    'Label Number: ', sprintf('%i',bnum)]};

% --- Executes when user attempts to close ParameterTool.
function ParameterTool_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to ParameterTool (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
delete(hObject);

