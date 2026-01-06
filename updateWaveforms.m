function flowCSF = updateWaveforms( ...
    flowCSF, ROITYPE, ...
    VplanesCSF, CSFSEG, ...
    INDEX, r, nframes)

% Extract and reshape ROI mask
roiCSF = CSFSEG.(ROITYPE)(INDEX, :);
imdim  = sqrt(numel(roiCSF));
roiCSF = reshape(roiCSF, imdim, imdim);

% Get CSF segment names (exclude tracking fields)
fns = setdiff(fieldnames(CSFSEG), ...
    {'coregTrack','segTrack','mthrTrack','sthrTrack','mansegTrack'}, ...
    'stable');

% Loop over time frames
for n = 1:nframes

    % --- Extract velocity planes ---
    for fld = ["x","y","z"]
        slice = squeeze(VplanesCSF.(fld)(INDEX,:,n));
        slice = reshape(slice, 2*r+1, []);
        v.(fld) = imresize(slice, [imdim imdim]); %#ok<AGROW>
    end

    % --- Total velocity ---
    vtotal = 0.1 * (v.x + v.y + v.z);

    % --- Apply mask ---
    maskedV = vtotal(roiCSF > 0);

    % --- Flow stats ---
    for i = 1:numel(fns)
        fn = fns{i};
        flowCSF.(fn).median(INDEX,n) = median(maskedV,'omitnan');
        flowCSF.(fn).mean(INDEX,n)   = mean(maskedV,'omitnan');
    end
end

end