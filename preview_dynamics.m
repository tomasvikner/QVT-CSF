function preview_dynamics(vxf, vyf, vzf, cxf, cyf, czf)
%PREVIEW_DYNAMICS  Blocking preview of central-slice CBF/CSF dynamics.
%   Top row: CBF (vx, vy, vz). Bottom row: CSF (vx, vy, vz).
%   Use slider / mouse wheel / left-right arrows to browse cardiac frames.
%   Execution pauses until this figure is closed.

tag = 'QVT_DynamicsPreview';
old = findall(0, 'Type', 'figure', 'Tag', tag);
if ~isempty(old)
    delete(old);
end

if ndims(vxf) < 4
    error('preview_dynamics:InputDims', 'Expected 4-D velocity arrays.');
end

nT = min([size(vxf,4), size(vyf,4), size(vzf,4), size(cxf,4), size(cyf,4), size(czf,4)]);
if nT < 1
    error('preview_dynamics:NoFrames', 'No temporal frames in velocity arrays.');
end

sz = size(vxf);
zc = max(1, round(sz(3) / 2));
t0 = 1;

f = figure('IntegerHandle', 'off', ...
    'HandleVisibility', 'callback', ...
    'NumberTitle', 'off', ...
    'WindowStyle', 'modal', ...
    'Name', sprintf('QVT: Dynamics preview (central z=%d/%d) — close to continue', zc, sz(3)), ...
    'Tag', tag, ...
    'Color', [0.08 0.08 0.08], ...
    'Position', [120 90 1120 760]);

t = tiledlayout(f, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact', ...
    'Units', 'normalized', 'Position', [0.04 0.18 0.92 0.78]);

ax = gobjects(6, 1);
for k = 1:6
    ax(k) = nexttile(t, k);
end

names = {'CBF Vx', 'CBF Vy', 'CBF Vz', 'CSF Vx', 'CSF Vy', 'CSF Vz'};

mx = @(A) max(abs(double(A(:))));
clims = [mx(vxf), mx(vyf), mx(vzf), mx(cxf), mx(cyf), mx(czf)];
clims(clims <= 0 | ~isfinite(clims)) = 1;

hIm = gobjects(6, 1);
for k = 1:6
    hIm(k) = imagesc(ax(k), zeros(sz(1), sz(2)));
    axis(ax(k), 'image');
    axis(ax(k), 'off');
    colormap(ax(k), 'parula');
    caxis(ax(k), [-clims(k), clims(k)]);
    title(ax(k), names{k}, 'Color', 'w', 'FontSize', 11);
end

lbl = uicontrol('Parent', f, 'Style', 'text', 'Units', 'normalized', ...
    'Position', [0.05 0.12 0.90 0.035], 'String', '', ...
    'BackgroundColor', get(f, 'Color'), 'ForegroundColor', [1 1 1], ...
    'HorizontalAlignment', 'center', 'FontSize', 10);

sld = uicontrol('Parent', f, 'Style', 'slider', 'Units', 'normalized', ...
    'Position', [0.12 0.075 0.76 0.035], 'Min', 1, 'Max', nT, 'Value', t0, ...
    'SliderStep', [1 / max(nT - 1, 1), 10 / max(nT - 1, 1)]);

uicontrol('Parent', f, 'Style', 'text', 'Units', 'normalized', ...
    'Position', [0.12 0.04 0.76 0.028], ...
    'String', 'Frame index (mouse wheel / left-right arrows also work)', ...
    'BackgroundColor', get(f, 'Color'), 'ForegroundColor', [1 1 1], ...
    'HorizontalAlignment', 'center');

set(f, 'CloseRequestFcn', @onClose, ...
    'WindowScrollWheelFcn', @onWheel, ...
    'KeyPressFcn', @onKey);
set(sld, 'Callback', @onSlider);

drawFrame(t0);
uiwait(f);

if isgraphics(f)
    set(f, 'CloseRequestFcn', []);
    delete(f);
end

    function drawFrame(ti)
        ti = max(1, min(nT, round(ti)));
        set(hIm(1), 'CData', vxf(:,:,zc,ti));
        set(hIm(2), 'CData', vyf(:,:,zc,ti));
        set(hIm(3), 'CData', vzf(:,:,zc,ti));
        set(hIm(4), 'CData', cxf(:,:,zc,ti));
        set(hIm(5), 'CData', cyf(:,:,zc,ti));
        set(hIm(6), 'CData', czf(:,:,zc,ti));
        set(lbl, 'String', sprintf('Central slice z=%d/%d | frame %d/%d', zc, sz(3), ti, nT));
        set(sld, 'Value', ti);
        drawnow limitrate;
    end

    function onSlider(src, ~)
        drawFrame(get(src, 'Value'));
    end

    function onWheel(~, ev)
        ti = round(get(sld, 'Value')) - sign(ev.VerticalScrollCount);
        drawFrame(ti);
    end

    function onKey(~, ev)
        ti = round(get(sld, 'Value'));
        switch ev.Key
            case {'rightarrow', 'uparrow'}
                ti = ti + 1;
            case {'leftarrow', 'downarrow'}
                ti = ti - 1;
            case 'home'
                ti = 1;
            case 'end'
                ti = nT;
            otherwise
                return;
        end
        drawFrame(ti);
    end

    function onClose(src, ~)
        uiresume(src);
    end

end
