function sol = solve_rect_levy_family(nu, k, n, boundaryKey)
%SOLVE_RECT_LEVY_FAMILY Analytic Levy-family modes for SS?? rectangular plates.
% x = 0,a are simply supported; y = 0,b take S/C/F according to boundaryKey.

a = 2.0;
b = 2.0;
Nplot = max(1201, 6*(n-2));
x = linspace(-1, 1, Nplot);
xp = x + 1.0;
y = x;
yp = y + 1.0;

mMax0 = max(18, ceil(sqrt(2*k)) + 12);
lambdaMax0 = max(400, 40*k);
Nscan = 5000;
rootAbsTol = 1e-8;
mergeTol = 1e-7;

modes = [];
coeffs = {};
for attempt = 1:3
    mMax = mMax0 + 6*(attempt-1);
    lambdaMax = lambdaMax0 * (1.0 + 0.75*(attempt-1));
    modes = [];
    coeffs = {};

    for m = 1:mMax
        alpha = m*pi/a;
        lower = alpha^2 + 1e-7;
        if lower >= lambdaMax
            continue;
        end
        lambdaGrid = linspace(lower, lambdaMax, Nscan);
        detVals = nan(size(lambdaGrid));
        for ii = 1:numel(lambdaGrid)
            detVals(ii) = levy_det(lambdaGrid(ii), m, nu, a, b, boundaryKey);
        end

        roots_m = [];
        for ii = 1:numel(lambdaGrid)-1
            l1 = lambdaGrid(ii);
            l2 = lambdaGrid(ii+1);
            f1 = detVals(ii);
            f2 = detVals(ii+1);
            if ~isfinite(f1) || ~isfinite(f2)
                continue;
            end
            useBracket = (sign(f1) ~= sign(f2));
            if ~useBracket && abs(f1) > 1e-6 && abs(f2) > 1e-6
                continue;
            end
            try
                r = fzero(@(lam) levy_det(lam, m, nu, a, b, boundaryKey), [l1, l2]);
                if isfinite(r) && r > lower && r <= lambdaMax
                    if abs(levy_det(r, m, nu, a, b, boundaryKey)) < rootAbsTol
                        roots_m(end+1) = r; %#ok<AGROW>
                    end
                end
            catch
            end
        end

        if isempty(roots_m)
            continue;
        end
        roots_m = sort(roots_m(:));
        roots_m = unique(round(roots_m, 8));
        roots_m = merge_close(roots_m, mergeTol);

        for ir = 1:numel(roots_m)
            lam = roots_m(ir);
            K = levy_matrix(lam, m, nu, a, b, boundaryKey);
            c = null_vector_from_matrix(K);
            modes = [modes; lam, m, ir]; %#ok<AGROW>
            coeffs{end+1} = c; %#ok<AGROW>
        end
    end

    if size(modes,1) >= k
        break;
    end
end

if isempty(modes)
    error('No Levy roots found for %s. Increase lambdaMax or mMax.', upper(boundaryKey));
end

[~, ord] = sort(modes(:,1), 'ascend');
modes = modes(ord, :);
coeffs = coeffs(ord);

kUse = min(k, size(modes,1));
modesU = cell(1, kUse);
lamDisp = modes(1:kUse, 1).';

for j = 1:kUse
    lam = modes(j,1);
    m = modes(j,2);
    c = coeffs{j};
    alpha = m*pi/a;
    [p, q] = pq_from_lambda(lam, alpha);
    Yv = levy_shape_vector(yp, p, q) * c;
    Xv = sin(m*pi*xp/a);
    U = real(Yv(:) * Xv(:).');
    U = canonicalize_mode(U);
    U = enforce_zero_edges(U, boundaryKey);
    utol = 1e-12 * max(1, max(abs(U(:)), [], 'omitnan'));
    U(abs(U) < utol) = 0;
    modesU{j} = U;
end

sol = struct('x', x, 'modesU', {modesU}, 'lamDisp', lamDisp);
end

function val = levy_det(lambda, m, nu, a, b, boundaryKey)
if lambda <= (m*pi/a)^2
    val = NaN;
    return;
end
K = levy_matrix(lambda, m, nu, a, b, boundaryKey);
rowNorms = vecnorm(K, 2, 2);
rowNorms(rowNorms < 1) = 1;
Ks = K ./ rowNorms;
val = det(Ks);
end

function K = levy_matrix(lambda, m, nu, a, b, boundaryKey)
alpha = m*pi/a;
[p, q] = pq_from_lambda(lambda, alpha);
rowY0 = boundary_rows(boundary_char(boundaryKey, 1), 0.0, p, q, alpha, nu);
rowYb = boundary_rows(boundary_char(boundaryKey, 2), b, p, q, alpha, nu);
K = [rowY0; rowYb];
end

function rows = boundary_rows(kind, y0, p, q, alpha, nu)
switch kind
    case 'S'
        rows = [RW(y0, p, q); RM(y0, p, q, alpha, nu)];
    case 'C'
        rows = [RW(y0, p, q); Rtheta(y0, p, q)];
    case 'F'
        rows = [RM(y0, p, q, alpha, nu); RV(y0, p, q, alpha, nu)];
    otherwise
        error('Unknown Levy boundary row kind: %s', kind);
end
end

function kind = boundary_char(boundaryKey, pos)
switch boundaryKey
    case 'ssss'
        chars = 'SS';
    case 'sscc'
        chars = 'CC';
    case 'ssff'
        chars = 'FF';
    case 'sssc'
        chars = 'SC';
    case 'sssf'
        chars = 'SF';
    case 'sscf'
        chars = 'CF';
    otherwise
        error('Unknown Levy boundary key: %s', boundaryKey);
end
kind = chars(pos);
end

function [p, q] = pq_from_lambda(lambda, alpha)
p = sqrt(alpha^2 + lambda);
q = sqrt(lambda - alpha^2);
end

function rows = levy_shape_vector(y, p, q)
rows = [cosh(p*y(:)), sinh(p*y(:)), cos(q*y(:)), sin(q*y(:))];
end

function out = RW(y, p, q)
out = [cosh(p*y), sinh(p*y), cos(q*y), sin(q*y)];
end

function out = Rtheta(y, p, q)
out = [p*sinh(p*y), p*cosh(p*y), -q*sin(q*y), q*cos(q*y)];
end

function out = RM(y, p, q, alpha, nu)
out = [(p^2 - nu*alpha^2)*cosh(p*y), ...
       (p^2 - nu*alpha^2)*sinh(p*y), ...
       -(q^2 + nu*alpha^2)*cos(q*y), ...
       -(q^2 + nu*alpha^2)*sin(q*y)];
end

function out = RV(y, p, q, alpha, nu)
cp = p*(p^2 - (2-nu)*alpha^2);
cq = q*(q^2 + (2-nu)*alpha^2);
out = [cp*sinh(p*y), cp*cosh(p*y), cq*sin(q*y), -cq*cos(q*y)];
end

function c = null_vector_from_matrix(K)
[~, ~, V] = svd(K, 'econ');
c = real(V(:,end));
if norm(c) < eps
    c = ones(4,1);
else
    c = c / norm(c);
end
if abs(c(1)) < 1e-12
    [~, idx] = max(abs(c));
else
    idx = 1;
end
if c(idx) < 0
    c = -c;
end
end

function v = merge_close(x, tol)
if isempty(x)
    v = x;
    return;
end
v = x(1);
for i = 2:numel(x)
    if abs(x(i) - v(end)) > tol
        v(end+1,1) = x(i); %#ok<AGROW>
    end
end
end

function U = enforce_zero_edges(U, boundaryKey)
% left/right correspond x=0,a and are always simply supported in Levy family
U(:,1) = 0;
U(:,end) = 0;
switch boundaryKey
    case {'ssss', 'sscc', 'sssc', 'sssf', 'sscf'}
        U(1,:) = 0;
end
switch boundaryKey
    case {'ssss', 'sscc', 'sssc'}
        U(end,:) = 0;
end
end

function v = canonicalize_mode(v)
[~, idx] = max(abs(v(:)));
if isempty(idx) || abs(v(idx)) < eps
    return;
end
if v(idx) < 0
    v = -v;
end
end
