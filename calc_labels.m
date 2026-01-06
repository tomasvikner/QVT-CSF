function [branchLabeled, Ntxt, txt] = calc_labels(branchLabeled, bnum, branchList, Labeltxt, value, average)

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

end