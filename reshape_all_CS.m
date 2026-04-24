function [CC, imdim] = reshape_all_CS(CSFSEG, CBFSEG, pindex)
CC = [];
fns = fieldnames(CSFSEG);
imdim = sqrt(size(CBFSEG,2)); % side length of cross-section
for i = 1:numel(fns) 
    fn = fns{i};
    if ~contains(fn, 'Track') % Don't update segmentation and coreg track
        CC.(fn) = reshape(double(CSFSEG.(fn)(pindex,:)),imdim,imdim);
    end
end
if sum(CC.madj(:)==0), CC.madj = CC.auto; end
end