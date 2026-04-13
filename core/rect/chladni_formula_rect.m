function chladni_formula_rect(nu, k, n, normalizeForDisplay, outputFolder, boundary)
%CHLADNI_FORMULA_RECT Rectangular / square plate modes.
%
% Preserved legacy branches:
%   FFFF -> original sparse free-edge formulation
%   CCCC -> original clamped finite-difference formulation
%   SSSS -> Navier exact modes
%
% Added analytic Levy-family branches:
%   SSCC, SSFF, SSSC, SSSF, SSCF

if nargin < 6 || isempty(boundary), boundary = 'free'; end
if nargin < 5 || isempty(outputFolder)
    outputFolder = fullfile(fileparts(mfilename('fullpath')), 'chladni_figures_output');
end
if ~exist(outputFolder, 'dir'), mkdir(outputFolder); end

meta = rect_boundary_meta(boundary);
nuStr = sprintf('%.6g', nu);

switch meta.solver_key
    case 'ffff'
        sol = solve_rect_free_sparse(nu, k, n);
    case 'ssss'
        sol = solve_rect_navier_ssss(k, n);
    case 'cccc'
        sol = solve_rect_clamped_fd_highres(k, n);
    case {'sscc', 'ssff', 'sssc', 'sssf', 'sscf'}
        sol = solve_rect_levy_family(nu, k, n, meta.solver_key);
    otherwise
        error('Unknown rectangular boundary condition: %s', boundary);
end

x = sol.x;
modesU = sol.modesU;
modesLamDisp = sol.lamDisp;

for j = 1:numel(modesU)
    U = modesU{j};

    fig = figure('Visible','off', 'Color', [1 1 1]);
    set(fig, 'InvertHardcopy', 'off');
    ax = gca;
    set(ax, 'Color', [1 1 1]);

    [Uf, climVal] = signed_field_for_display(U, normalizeForDisplay);
    imagesc(ax, x, x, Uf);
    set(ax, 'YDir', 'normal');
    axis(ax, 'equal');
    axis(ax, [-1 1 -1 1]);

    color_bar(ax, 'Location','eastoutside', 'Interpreter','latex', 'Limits',[-climVal climVal]);
    hold(ax, 'on');
    draw_rect_nodal_lines(ax, x, U, meta.solver_key);

    apply_latex_formatting(fig, ax);
    set(ax, 'Layer', 'top');
    title(ax, sprintf('$\\nu=%.6g,\\ %s\\ (\\mathrm{mode}=%d),\\ \\Lambda=%.4g$', ...
        nu, meta.title_tag, j, modesLamDisp(j)), 'Interpreter','latex');

    filename = fullfile(outputFolder, sprintf('rect-%s-%s-%02d.png', meta.file_tag, nuStr, j));
    print(fig, '-dpng', filename);
    close(fig);
    fprintf('  Graphics saved as "%s"\n', filename);
end
end

function draw_rect_nodal_lines(ax, x, U, boundaryKey)
umax = max(abs(U(:)), [], 'omitnan');
if ~isfinite(umax) || umax < eps
    umax = 1.0;
end
U2 = U;
U2(abs(U2) < 1e-12 * umax) = 0;
contour(ax, x, x, U2, [0 0], 'k-', 'LineWidth', 1.0);

[leftZero, rightZero, bottomZero, topZero] = boundary_zero_edges(boundaryKey);
if leftZero
    plot(ax, [-1 -1], [-1 1], 'k-', 'LineWidth', 1.0);
end
if rightZero
    plot(ax, [1 1], [-1 1], 'k-', 'LineWidth', 1.0);
end
if bottomZero
    plot(ax, [-1 1], [-1 -1], 'k-', 'LineWidth', 1.0);
end
if topZero
    plot(ax, [-1 1], [1 1], 'k-', 'LineWidth', 1.0);
end
end

function [leftZero, rightZero, bottomZero, topZero] = boundary_zero_edges(boundaryKey)
leftZero = false;
rightZero = false;
bottomZero = false;
topZero = false;

switch boundaryKey
    case 'ffff'
        return;
    case 'cccc'
        leftZero = true; rightZero = true; bottomZero = true; topZero = true;
    case 'ssss'
        leftZero = true; rightZero = true; bottomZero = true; topZero = true;
    case 'sscc'
        leftZero = true; rightZero = true; bottomZero = true; topZero = true;
    case 'ssff'
        leftZero = true; rightZero = true;
    case 'sssc'
        leftZero = true; rightZero = true; bottomZero = true; topZero = true;
    case 'sssf'
        leftZero = true; rightZero = true; bottomZero = true;
    case 'sscf'
        leftZero = true; rightZero = true; bottomZero = true;
end
end
