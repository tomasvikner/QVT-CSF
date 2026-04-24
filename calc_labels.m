function [MD, txt] = calc_labels(MD, bnum, branchList, value, average)

% Put the number labels on the CenterlinePlot if new branch
if MD.branchLabeled ~= bnum
    delete(MD.Ntxt)
    MD.branchLabeled = bnum;
    index_branch = branchList(:,4) == MD.branchLabeled;
    branchActual = branchList(index_branch,1:3);
    textint = 0:5:length(branchActual)-1;
    numString_val = num2str(textint);
    numString_val = strsplit(numString_val);
    MD.Ntxt = text(branchActual(textint+1,1),branchActual(textint+1,2), ...
        branchActual(textint+1,3),numString_val,'Color','w','FontSize',10,...
        'HitTest','off','PickableParts','none','Parent',fig.CurrentAxes);
end

% Get branch indices and current label point
branchActual = branchList(branchList(:,4) == MD.branchLabeled,5);
CurrentNum = find(branchList(pindex,5)==branchActual)-1;

% Update cursor text
Labeltxt = MD.Labeltxt;
ptLabel = MD.PointLabel;
txt = {['Point Label: ' , ptLabel , sprintf('\n'), ...
    Labeltxt{1,1}, sprintf('%0.3f',value),Labeltxt{1,2}, sprintf('\n'), ...
    Labeltxt{2,1},sprintf('%0.3f',mean(average)),Labeltxt{2,2},sprintf('\n'), ...
    'Current Branch #: ',sprintf('%i',CurrentNum),sprintf('\n') ...
    'Label Number: ', sprintf('%i',bnum)]};

end