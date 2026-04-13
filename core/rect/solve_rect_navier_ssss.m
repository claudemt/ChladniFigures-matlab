function sol = solve_rect_navier_ssss(k, n)
%SOLVE_RECT_NAVIER_SSSS Exact Navier modes for the square plate.

a = 2.0;
b = 2.0;
Nplot = max(1201, 6*(n-2));
x = linspace(-1, 1, Nplot);
xp = x + 1.0;
yp = xp;

mMax = max(12, ceil(sqrt(k)) + 8);
pairs = zeros(mMax*mMax, 2);
lam = zeros(mMax*mMax, 1);
t = 0;
for m = 1:mMax
    for nn = 1:mMax
        t = t + 1;
        pairs(t,:) = [m nn];
        lam(t) = (m*pi/a)^2 + (nn*pi/b)^2;
    end
end
pairs = pairs(1:t,:);
lam = lam(1:t);
[lam, ord] = sort(lam, 'ascend');
pairs = pairs(ord,:);

kUse = min(k, numel(lam));
modesU = cell(1, kUse);
lamDisp = lam(1:kUse).';
lamDisp(lamDisp < 1e-12) = 0;

for j = 1:kUse
    m = pairs(j,1);
    nn = pairs(j,2);
    X = sin(m*pi*xp/a);
    Y = sin(nn*pi*yp/b);
    U = Y(:) * X(:).';
    U(1,:) = 0; U(end,:) = 0; U(:,1) = 0; U(:,end) = 0;
    modesU{j} = U;
end

sol = struct('x', x, 'modesU', {modesU}, 'lamDisp', lamDisp);
end
