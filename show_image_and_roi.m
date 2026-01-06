function show_img_and_roi(ax, img, ttl, mask, clim)

if nargin < 5 || isempty(clim)
    imshow(img, [], 'InitialMagnification', 'fit', 'Parent', ax);
else
    imshow(img, clim, 'InitialMagnification', 'fit', 'Parent', ax);
end

title(ax, ttl, 'FontSize', 13);

if ~isempty(mask)
    visboundaries(ax, mask, 'LineWidth', 1);
end

end