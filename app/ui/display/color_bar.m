function out = color_bar(ax, varargin)
%COLOR_BAR Visible-spectrum colormap + unified colorbar rules.

p = inputParser;
p.addRequired('ax', @(h) ishghandle(h,'axes'));
p.addParameter('N', 256, @(v) isnumeric(v) && isscalar(v) && v>=2);
p.addParameter('UseVisibleSpectrum', true, @(v) islogical(v) || (isnumeric(v)&&isscalar(v)));
p.addParameter('Colormap', [], @(v) isempty(v) || (isnumeric(v)&&size(v,2)==3));
p.addParameter('Limits', [], @(v) isempty(v) || (isnumeric(v)&&numel(v)==2 && all(isfinite(v))));
p.addParameter('AutoSymmetric', false, @(v) islogical(v) || (isnumeric(v)&&isscalar(v)));
p.addParameter('Data', [], @(v) isempty(v) || isnumeric(v));
p.addParameter('NormalizeToUnit', false, @(v) islogical(v) || (isnumeric(v)&&isscalar(v)));
p.addParameter('CreateColorbar', true, @(v) islogical(v) || (isnumeric(v)&&isscalar(v)));
p.addParameter('Location', 'eastoutside', @(s) ischar(s) || isstring(s));
p.addParameter('Interpreter', 'latex', @(s) ischar(s) || isstring(s));
p.addParameter('Label', '', @(s) ischar(s) || isstring(s));
p.addParameter('LabelInterpreter', '', @(s) ischar(s) || isstring(s));
p.addParameter('Ticks', [], @(v) isempty(v) || isnumeric(v));
p.addParameter('TickLabels', [], @(v) isempty(v) || isnumeric(v) || isstring(v) || iscellstr(v));
p.parse(ax, varargin{:});
opt = p.Results;

if ~isempty(opt.Colormap)
    cmap = opt.Colormap;
else
    if opt.UseVisibleSpectrum
        cmap = viscolormap_local(opt.N);
    else
        cmap = parula(opt.N);
    end
end
colormap(ax, cmap);

climApplied = [];
if ~isempty(opt.Limits)
    clim(ax, opt.Limits);
    climApplied = opt.Limits;
elseif opt.AutoSymmetric
    data = opt.Data;
    if isempty(data)
        data = try_get_axes_cdata(ax);
    end
    if ~isempty(data)
        if opt.NormalizeToUnit
            s = max(abs(data(:)), [], 'omitnan');
            if ~isfinite(s) || s < eps, s = 1; end
            data = data ./ s;
        end
        m = max(abs(data(:)), [], 'omitnan');
        if ~isfinite(m) || m < eps, m = 1; end
        clim(ax, [-m m]);
        climApplied = [-m m];
    end
end

cb = [];
if opt.CreateColorbar
    cb = colorbar(ax, char(opt.Location));
    set(cb, 'TickLabelInterpreter', char(opt.Interpreter));
    if ~isempty(opt.Ticks), cb.Ticks = opt.Ticks; end
    if ~isempty(opt.TickLabels), cb.TickLabels = opt.TickLabels; end
    if strlength(string(opt.Label)) > 0
        cb.Label.String = char(opt.Label);
        if strlength(string(opt.LabelInterpreter)) == 0
            cb.Label.Interpreter = char(opt.Interpreter);
        else
            cb.Label.Interpreter = char(opt.LabelInterpreter);
        end
    end
end

out = struct('cb', cb, 'cmap', cmap, 'clim', climApplied);
end

function data = try_get_axes_cdata(ax)
data = [];
kids = ax.Children;
for i = 1:numel(kids)
    h = kids(i);
    if isprop(h, 'CData')
        cd = h.CData;
        if isnumeric(cd) && ~isempty(cd)
            data = cd;
            return;
        end
    end
    if isprop(h, 'ZData')
        zd = h.ZData;
        if isnumeric(zd) && ~isempty(zd)
            data = zd;
            return;
        end
    end
end
end

function cmapOut = viscolormap_local(N)
if nargin < 1, N = 256; end
lambda = linspace(380, 780, N);
rgb = zeros(N, 3);
for ii = 1:N
    l = lambda(ii);
    if l >= 380 && l < 440
        r = -(l - 440) / 60; g = 0; b = 1;
    elseif l >= 440 && l < 490
        r = 0; g = (l - 440) / 50; b = 1;
    elseif l >= 490 && l < 510
        r = 0; g = 1; b = -(l - 510) / 20;
    elseif l >= 510 && l < 580
        r = (l - 510) / 70; g = 1; b = 0;
    elseif l >= 580 && l < 645
        r = 1; g = -(l - 645) / 65; b = 0;
    elseif l >= 645 && l <= 780
        r = 1; g = 0; b = 0;
    else
        r = 0; g = 0; b = 0;
    end

    if l >= 380 && l < 420
        f = 0.3 + 0.7*(l - 380)/40;
    elseif l >= 420 && l <= 700
        f = 1.0;
    elseif l > 700 && l <= 780
        f = 0.3 + 0.7*(780 - l)/80;
    else
        f = 0.0;
    end

    gamma = 0.8;
    rgb(ii, :) = (f .* [r g b]) .^ gamma;
end
cmapOut = max(min(rgb, 1), 0);
end
