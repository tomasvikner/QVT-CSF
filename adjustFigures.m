%% Edit figures before opening GUI  
addpath QVT; 
close all; 

% Load the .fig file
openfig('uu/paramMap.fig', 'visible');
fig = gcf;

axesHandles = findall(fig, 'Type', 'axes');
numSubplots = numel(axesHandles);
% Resize the figure
movegui(fig, 'east');
fig.Position(2) = fig.Position(2) * 0.15;
fig.Position(4) = fig.Position(4) * 2.0;
fig.Position(3) = fig.Position(3) * 1.0; 
fig.Position(1) = fig.Position(1) -0.10 * fig.Position(1);

csahi = [1 3 4 5];
ahlist = zeros(5, 4);  

for i = csahi
    ah = axesHandles(i);
    ah.Position(4) = ah.Position(4) - 0.5 * ah.Position(4); 
    ahlist(i, :) = ah.Position; 
end

% Add new axes with specified positions
newNames = {'CSF-Mag', 'COB', 'CXC', 'CSF-CD'};
newAxesHandles = gobjects(4, 1); % Pre-allocate for new axes handles
for i = 1:4
    j = csahi(i);
    newpos = ahlist(j, :); newpos(2) = newpos(2) + newpos(4);
    % newname = ['New Axis ' num2str(i)]; 
    newname = newNames{5-i};
    newAxesHandles(i) = axes('Parent', fig, 'Position', newpos);
    title(newAxesHandles(i), newname);
end

% newname = 'csfwaveform';
% newAxesHandles(5) = axes('Parent', fig, 'Position', newpos);

% Add some manual adjustments 
fig = gcf;
axesHandles = findall(fig, 'Type', 'axes');
axesHandles(6).Position(2) = axesHandles(6).Position(2) - 0.05; % Shift down waveform plot 
axesHandles(6).Position(4) = axesHandles(6).Position(4) - 0.05; % Shift down waveform plot 
axesHandles(6).Position(3) = axesHandles(6).Position(3) - 0.06; % Add space for CSF y-label

% Shift down 4 cross sections 
for i = 1:4
    j = csahi(i) + 4; 
    % [i, j]

    ah = axesHandles(j);
    ah.Position(2) = ah.Position(2) - 0.08; 
    ah.Position(4) = ah.Position(4) + 0.03;

    ah = axesHandles(i);
    ah.Position(2) = ah.Position(2) - 0.03; 
    ah.Position(4) = ah.Position(4) + 0.03;
    ah.XTick = [];
    ah.YTick = [];
    
end

newfig = fig; 
allObjects = findall(newfig);

% Additional code to find the UIControl with 'parameter_choice' tag
foundParamChoice = false;

% Look for UIControl named 'parameter_choice' and append new string
for obj = allObjects'
    if strcmp(get(obj, 'Type'), 'uicontrol') && strcmp(get(obj, 'Tag'), 'parameter_choice')
        foundParamChoice = true;
        
        % Retrieve the current strings in the 'parameter_choice' UIControl
        paramStrings = get(obj, 'String');
        
        % Check if paramStrings is a cell array of strings
        if iscell(paramStrings) && length(paramStrings) == 8
            % Append the new string to the cell array
            paramStrings{end+1} = 'CSF Coupling';
            
            % Update the 'String' property with the modified list
            set(obj, 'String', paramStrings);
            disp('New string appended successfully!');
        else
            disp('Expected an 8-element cell array for the parameter choice strings.');
        end
        break;  % Exit loop once the target object is found
    end
end

% Display message if the parameter choice UIControl was not found
if ~foundParamChoice
    disp('UIControl with Tag "parameter_choice" not found.');
end

% Locate the Manual Segmentation button
manualSegButton = findall(fig, 'Style', 'pushbutton', 'String', 'Manual Segmentation');

if ~isempty(manualSegButton)
    % Reduce size and shift left
    manualSegButton.Position(3) = manualSegButton.Position(3) * 0.6;
    manualSegButton.Position(4) = manualSegButton.Position(4) * 1.5;
    manualSegButton.Position(2) = manualSegButton.Position(2) + 0.00;
    manualSegButton.Position(1) = manualSegButton.Position(1) - 0.025;
    
    % Add "Manual Coreg" button next to it
    manualCoregButton = uicontrol('Style', 'pushbutton', 'Parent', fig, ...
        'String', 'Manual Coreg', 'FontSize', 10, ...
        'Units', 'normalized', 'Position', ...
        [manualSegButton.Position(1) + manualSegButton.Position(3) + 0.02, ...
         manualSegButton.Position(2), ...
         manualSegButton.Position(3), manualSegButton.Position(4)], ...
        'Tag', 'manual_coreg', 'Callback', @(hObject, eventdata)paramMap('manualCoreg_Callback', hObject, eventdata, guidata(hObject)));
end

% Make a waveform toggling button

wfb2 = uicontrol('Style', 'pushbutton', 'Parent', fig, ...
        'String', 'PC-1 WF', 'FontSize', 10, ...
        'Units', 'normalized', 'Position', [0.88, 0.18 0.07 0.04], ...
        'Tag', 'manual_coreg', 'Callback', @(hObject, eventdata)paramMap('togglePCA_Callback', hObject, eventdata, guidata(hObject)));

wfb1 = uicontrol('Style', 'pushbutton', 'Parent', fig, ...
        'String', 'Vel STD WF', 'FontSize', 10, ...
        'Units', 'normalized', 'Position', [0.81, 0.18 0.07 0.04], ...
        'Tag', 'manual_coreg', 'Callback', @(hObject, eventdata)paramMap('toggleCnoB_Callback', hObject, eventdata, guidata(hObject)));

wfb3 = uicontrol('Style', 'pushbutton', 'Parent', fig, ...
        'String', 'CUBE WF', 'FontSize', 10, ...
        'Units', 'normalized', 'Position', [0.74, 0.18 0.07 0.04], ...
        'Tag', 'manual_coreg', 'Callback', @(hObject, eventdata)paramMap('toggleCUBE_Callback', hObject, eventdata, guidata(hObject)));

% Positioning of 8x segmentation boxes 
for i = [1, 2, 3, 4, 5, 7, 8, 9]
    ahi = axesHandles(i);
    ahi.Position(4) = 0.18;
    ahi.Position(3) = 0.22;
    [ahi.Position(1), ahi.Position(2)]
end
x1 = 0.04; 
x2 = 0.27;
x3 = 0.50; 
x4 = 0.73; 

y4 = 0.735; y3 = 0.53;
axesHandles(1).Position(1) = x2; 
axesHandles(1).Position(2) = y4; 
axesHandles(2).Position(1) = x3; 
axesHandles(2).Position(2) = y4; 
axesHandles(3).Position(1) = x4; 
axesHandles(3).Position(2) = y4; 
axesHandles(4).Position(1) = x1; 
axesHandles(4).Position(2) = y4; 
axesHandles(5).Position(1) = x1;  
axesHandles(5).Position(2) = y3; 
axesHandles(7).Position(1) = x4; 
axesHandles(7).Position(2) = y3; 
axesHandles(8).Position(1) = x3; 
axesHandles(8).Position(2) = y3; 
axesHandles(9).Position(1) = x2; 
axesHandles(9).Position(2) = y3; 

newNames = {'CUBE Anti-FLAIR', 'CSF velocity STD', 'CSF vel. (time resolved)', 'CSF magnitude', ...
    'CBF magnitude', 'CBF vel. (time resolved)', 'Complex difference (CD)', 'T1-w. MP-RAGE'};
c = 0; 
for i = [1, 2, 3, 4, 5, 7, 8, 9]
    c = c + 1; 
    newname = newNames{c};
    ttl = title(axesHandles(i), newname);
    ttl.FontSize = 12; 
end

% New Pos Mars 25
for obj = allObjects'
    if isprop(obj, 'String')  % Check if object has a 'String' property
        disp(['Object Type: ', get(obj, 'Type')]);

        if contains(get(obj, 'String'), 'Magnitude')
            disp('Match found!');
            set(obj, 'String', sprintf('CSF-Mag (BICS)\nCBF-Mag (CBF)'));
            set(obj, 'String', sprintf(''));
        elseif contains(get(obj, 'String'), 'Complex Difference')
            disp('Match found!');
            set(obj, 'String', sprintf('T2 CUBE (BCSF)\nMPRAGE (CBF)'));
            set(obj, 'String', sprintf(''));
        elseif contains(get(obj, 'String'), 'Time-Averaged')
            disp('Match found!');
            set(obj, 'String', sprintf('CSF-SD (CNOB)\nCBF-CD (CBF)')); % SD is more of interest than TA velocity 
            set(obj, 'String', sprintf(''));
        elseif contains(get(obj, 'String'), 'Time-Resolved')
            disp('Match found!');
            set(obj, 'String', sprintf('Vel. Time-Resolved\nCSF (top); CBF (bot)'));
            set(obj, 'String', sprintf(''));
        end

    end
end

% Sliders 
areaSlider = allObjects(10);

% Get position of existing areaSlider
sliderPos = areaSlider.Position;

initialValues = 0.35;

% Optional: Save the resized figure
savefig(newfig, 'QVT/paramMap.fig');
disp('Saved new Figure'); 

