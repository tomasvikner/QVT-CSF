function hpatch = show_iso_seg(ax, segment, hpatch)

if ~isempty(hpatch) && isgraphics(hpatch)
    delete(hpatch);
end

axes(ax); hold(ax,'on');

hpatch = patch(isosurface(permute(segment,[2 1 3]),0.5), ...
               'FaceAlpha',0);

reducepatch(hpatch,0.7);

set(hpatch, ...
    'FaceColor','white', ...
    'EdgeColor','none', ...
    'PickableParts','none');

set(ax,'Color','black','ZDir','reverse');
axis(ax,'off','tight','vis3d');
daspect(ax,[1 1 1]);
view(ax,[0 0 1]);

camlight(ax,'headlight');
lighting(ax,'gouraud');

hold(ax,'off');

end