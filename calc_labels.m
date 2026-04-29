function [MD, txt] = calc_labels(MD, bnum, branchList, value, average, pindex)
% pindex: centerline point index (required for branch label text)

global fig %#ok<GVMIS>

% Put the number labels on the CenterlinePlot if new branch
if MD.branchLabeled ~= bnum
    delete(MD.Ntxt)
    MD.branchLabeled = bnum;
    index_branch = branchList(:,4) == MD.branchLabeled;
    branchActual = branchList(index_branch,1:3);
    textint = 0:5:length(branchActual)-1;
    numString_val = num2str(textint);
    numString_val = strsplit(numString_val);
    if ~isempty(fig) && isgraphics(fig) && isgraphics(fig.CurrentAxes)
        parentAx = fig.CurrentAxes;
    else
        parentAx = [];
    end
    if isempty(parentAx)
        parentAx = gca;
    end
    MD.Ntxt = text(branchActual(textint+1,1),branchActual(textint+1,2), ...
        branchActual(textint+1,3),numString_val,'Color','w','FontSize',10,...
        'HitTest','off','PickableParts','none','Parent',parentAx);
end

% Get branch indices and current label point
if size(branchList, 2) >= 5
    branchActual = branchList(branchList(:,4) == MD.branchLabeled, 5);
    ix = find(branchList(pindex, 5) == branchActual, 1);
    if isempty(ix)
        CurrentNum = NaN;
    else
        CurrentNum = ix - 1;
    end
else
    CurrentNum = NaN;
end

if isnan(CurrentNum)
    curBranchStr = 'n/a';
else
    curBranchStr = sprintf('%i', CurrentNum);
end

% Update cursor text
Labeltxt = MD.Labeltxt;
ptLabel = MD.PointLabel;
txt = {['Point Label: ' , ptLabel , sprintf('\n'), ...
    Labeltxt{1,1}, sprintf('%0.3f',value),Labeltxt{1,2}, sprintf('\n'), ...
    Labeltxt{2,1},sprintf('%0.3f',mean(average)),Labeltxt{2,2},sprintf('\n'), ...
    'Current Branch #: ', curBranchStr, sprintf('\n') ...
    'Label Number: ', sprintf('%i',bnum)]};

end