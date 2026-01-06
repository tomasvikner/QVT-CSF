function hscatter = area_thresh(handles, hscatter, CLVALS, branchList, LogPoints)
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