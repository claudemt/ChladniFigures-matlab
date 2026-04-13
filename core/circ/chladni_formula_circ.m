function chladni_formula_circ(nu, k, n, normalizeForDisplay, outputFolder, boundary)
%CHLADNI_FORMULA_CIRC Circular plate (Bessel / modified Bessel roots).
% boundary: 'clamped' | 'simply' | 'free'
% Output: circ-F/S/C-<nu>-<mode>.png

if nargin < 6 || isempty(boundary), boundary = 'free'; end
if nargin < 5 || isempty(outputFolder)
    outputFolder = fullfile(fileparts(mfilename('fullpath')), 'chladni_figures_output');
end
if ~exist(outputFolder, 'dir'), mkdir(outputFolder); end

boundary = char(lower(string(boundary)));
nuStr = sprintf('%.6g', nu);

a = 1.0;
x = linspace(-a, a, n);
y = x;
[X, Y] = meshgrid(x, y);
Rr = hypot(X, Y);
TH = atan2(Y, X);
mask = (Rr <= a);

mMax = max(12, ceil(sqrt(2*k)) + 8);
betaMax = max(80, 40 + 10*sqrt(k));
step = 0.5;

modes = []; % rows: [beta, m, n_idx, C]
for m = 0:mMax
    roots_m = [];
    b0 = 1e-6;
    f0 = feval_boundary_char(b0, m, nu, boundary);
    for b1 = b0+step:step:betaMax
        f1 = feval_boundary_char(b1, m, nu, boundary);
        if isfinite(f0) && isfinite(f1) && sign(f0) ~= sign(f1)
            try
                br = fzero(@(bb) feval_boundary_char(bb, m, nu, boundary), [b0, b1]);
                if isfinite(br) && br > 1e-4
                    roots_m(end+1) = br; %#ok<AGROW>
                end
            catch
            end
        end
        b0 = b1;
        f0 = f1;
    end

    if isempty(roots_m), continue; end
    roots_m = unique(round(roots_m, 10));
    for ni = 1:numel(roots_m)
        beta = roots_m(ni);
        C = -besselj(m, beta) / besseli(m, beta);
        modes = [modes; beta, m, ni, C]; %#ok<AGROW>
    end
end

if isempty(modes)
    error('No circular eigen-roots found. Increase betaMax or mMax.');
end

[~, idxSort] = sort(modes(:,1).^2, 'ascend');
modes = modes(idxSort, :);
kUse = min(k, size(modes,1));
rr = Rr / a;
bcLetter = boundary_short(boundary);

for i = 1:kUse
    beta = modes(i,1);
    m = modes(i,2);
    ni = modes(i,3);
    C = modes(i,4);

    U = nan(n, n);
    core = besselj(m, beta*rr(mask)) + C * besseli(m, beta*rr(mask));
    if m > 0
        core = core .* cos(m * TH(mask));
    end
    U(mask) = core;

    fig = figure('Visible','off', 'Color', [1 1 1]);
    set(fig, 'InvertHardcopy', 'off');
    set(fig, 'Renderer', 'opengl');
    ax = gca; set(ax, 'Color', [1 1 1]);

    [Uf, climVal] = signed_field_for_display(U, normalizeForDisplay);
    UfImg = Uf;
    UfImg(~mask) = 0;
    hImg = imagesc(x, y, UfImg);
    set(ax, 'YDir', 'normal');
    axis(ax, 'equal');
    axis(ax, [-a a -a a]);
    set(hImg, 'AlphaData', double(mask), 'AlphaDataMapping', 'none');

    color_bar(ax, 'Location','eastoutside', 'Interpreter','latex', 'Limits',[-climVal climVal]);
    hold(ax, 'on');
    contour(x, y, U, [0 0], 'k-', 'LineWidth', 1.0);

    apply_latex_formatting(fig, ax);
    lam_i = beta^2 / a^2;
    title(ax, sprintf('$\\nu=%.6g,\\ %s\\ (m=%d,\\ n=%d),\\ \\Lambda=%.4g$', ...
        nu, bcLetter, m, ni, lam_i), 'Interpreter', 'latex');

    filename = fullfile(outputFolder, sprintf('circ-%s-%s-%02d.png', bcLetter, nuStr, i));
    print(fig, '-dpng', filename);
    close(fig);
    fprintf('  Graphics saved as "%s"\n', filename);
end
end

function letter = boundary_short(bc)
switch bc
    case 'clamped', letter = 'C';
    case 'simply',  letter = 'S';
    case 'free',    letter = 'F';
    otherwise, error('Unknown boundary condition: %s', bc);
end
end

function val = feval_boundary_char(beta, m, nu, bc)
switch bc
    case 'clamped'
        J = besselj(m, beta); I = besseli(m, beta);
        Jp1 = besselj(m+1, beta); Ip1 = besseli(m+1, beta);
        val = J * Ip1 + I * Jp1;
    case 'simply'
        J = besselj(m, beta); I = besseli(m, beta);
        Jp1 = besselj(m+1, beta); Ip1 = besseli(m+1, beta);
        val = Jp1*I + Ip1*J - (2*beta/(1-nu))*(J*I);
    case 'free'
        J = besselj(m, beta); I = besseli(m, beta);
        Jd = besselj_deriv(m, beta); Id = besseli_deriv(m, beta);
        r11 = beta^2*J + (1-nu)*(beta*Jd - m^2*J);
        r12 = -beta^2*I + (1-nu)*(beta*Id - m^2*I);
        r21 = beta^3*Jd + m^2*(1-nu)*(beta*Jd - J);
        r22 = -beta^3*Id + m^2*(1-nu)*(beta*Id - I);
        s1 = max(1, hypot(r11, r12));
        s2 = max(1, hypot(r21, r22));
        val = (r11/s1)*(r22/s2) - (r12/s1)*(r21/s2);
    otherwise
        error('Unknown boundary condition: %s', bc);
end
end

function Jd = besselj_deriv(m, x)
if m == 0
    Jm1 = -besselj(1, x);
    Jp1 =  besselj(1, x);
else
    Jm1 = besselj(m-1, x);
    Jp1 = besselj(m+1, x);
end
Jd = 0.5 * (Jm1 - Jp1);
end

function Id = besseli_deriv(m, x)
if m == 0
    Im1 = besseli(1, x);
    Ip1 = besseli(1, x);
else
    Im1 = besseli(m-1, x);
    Ip1 = besseli(m+1, x);
end
Id = 0.5 * (Im1 + Ip1);
end
