% 01/05/26: move this to stand alone function 
function hscatter = invert_area(hObject, hscatter, branchList, CLVALS, handles)
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