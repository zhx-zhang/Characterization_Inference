% run EM algorithm to find the solution
% this will produce most of the plots shown in the paper

% intGaussGrid.mat: pre-computed data for fast discrete Gaussian calculation
global M LAM MU VA
load data/intGaussGrid
VA = S;
LAM = L;

load data/inference_input

A_exp = cat(1, A{:});
y_exp = cat(1, y{:});
m0 = A_exp \ y_exp; % nominal fitting, naive approach 

n = length(y);
pn = length(m0);
AA = zeros(6, pn, n);
yy = cellfun(@(y, e) y./e.^2, y, e_y, 'UniformOutput', false);
dd = cellfun(@(d, e) d./e.^2, d_y, e_y, 'UniformOutput', false);
IS = cellfun(@pinv, S, 'UniformOutput', false);
LL = cellfun(@(e) diag(1./e.^2), e_y, 'UniformOutput', false);
li0 = cellfun(@(d, e) sum((d./e).^2), d_y, e_y);
ml0 = cellfun(@(d, y, e) sum(d.*y./e.^2), d_y, y, e_y);

for i = 1 : n
    AA(:, :, i) = IS{i} * A{i};
end

% some settings for the algorithm
N_ITER = 20000; % maximal iteration
tol = 1e-4; % convergence test tolerance
Sig0 = diag([27, 0.8, 4.5, 0.1, 0.01, 0.1].^2); % initial guess
m = m0; % intial guess
l_pr = .5; % a small prior/regularization for discrete Gaussian, helps with stability

% initialization
n_p = 4;
Np = n + n_p - 4;
IR = cell(n, 1);
ISR = cell(n, 1);
K = zeros(6, 6, n);
mu = zeros(6, 1);
V = inv(Sig0);
ni = zeros(n, 1);

% EM iteration
for iter = 1 : N_ITER
    % update q & li
    Lq = zeros(pn); Lqm = zeros(pn, 1); li = li0;
    for i = 1 : n
	ISi = IS{i};
	IRi = inv(ISi.'*V*ISi + LL{i});
	ISRi = ISi*IRi;
	Ki = ISRi*ISi.';
	Ai = AA(:, :, i);
	AVi = Ai.'*V;
	Lq = Lq + AVi*Ai;
	Lqm = Lqm + AVi * (Ki*AVi.'*m + ISRi*(yy{i} - ni(i)*dd{i}));
	li(i) = li(i) - dd{i}.'*IRi*dd{i};
	K(:, :, i) = Ki; IR{i} = IRi; ISR{i} = ISRi;
    end
    m = Lq \ Lqm;

    % update mi
    ml = ml0;
    for i = 1 : n
	ml(i) = ml(i) - dd{i}.'*(ISR{i}.'*V*AA(:, :, i)*m + IR{i}*yy{i});
    end
    % prior on mi & li
    li = li + l_pr;
    mi = ml./li;
    [ni, vi] = intGaussian(mi, li);

    % update V
    K_sum = sum(K, 3);
    IV = K_sum;% + n_p*Sig0; for prior version
    for i = 1 : n
	Ki = K(:, :, i);
	ISRi = ISR{i};
	mu(:, i) = Ki*V*AA(:, :, i)*m + ISRi*(yy{i} - ni(i)*dd{i});
	r = mu(:, i) - AA(:, :, i)*m;
	sd = ISRi*dd{i};
	IV  = IV + r*r.' + vi(i)*sd*sd.';
	K(:, :, i) = Ki + vi(i)*sd*sd.';
    end
    V1 = inv(IV/Np);

    % calculate the increment, determine convergence
    u0 = -1/2*sum(dot(IV, V)) + Np/2*log(det(V));
    u1 = -Np*6/2 + Np/2*log(det(V1));
    if u1 - u0 < tol
	break;
    else
	V = V1;
    end
end

% convert to the p-space for the posterior
dp = zeros(6, n);
dp0 = zeros(6, n);
n_ellip = 20;
[X, Y, Z] = sphere(n_ellip);
u = [X(:) Y(:) Z(:)]';

for i = 1 : n
    % plot the naive results
    dp0(:, i) = IS{i} * y{i} - AA(:, :, i)*m0;
    figure(1);
    if i ~= 29 % this is an outlier that's very wrong
        plot3(dp0(1, i), dp0(2, i), dp0(3, i), 'k.', 'MarkerSize', 12); hold on;
    end

    % bayesian results
    dp(:, i) = mu(:, i) - AA(:, :, i)*m;
    figure(2);
    plot3(dp(1, i), dp(2, i), dp(3, i), 'k.', 'MarkerSize', 12); hold on;
    % this is for plotting the distribution ellipse
    [COEFF, latent] = eig(K(1:3, 1:3, i));
    alpha = max(min(0.01/sqrt(prod(diag(latent))), 0.1), 0.02);
    w = 4 * (COEFF * sqrt(latent)) * u;
    z = repmat(dp(1:3, i), 1, (n_ellip+1)^2) + w;
    surf(reshape(z(1, :), n_ellip+1, n_ellip+1), ...
        reshape(z(2, :), n_ellip+1, n_ellip+1), ...
        reshape(z(3, :), n_ellip+1, n_ellip+1), ...
        'EdgeColor', 'none', 'FaceColor', '#0072BD', ...
        'FaceAlpha', alpha);
end

for i = 1 : 2
    figure(i);
    hold off;
    xlabel('\Deltaw (nm)'); ylabel('\DeltaT (nm)'); zlabel('\Deltah (nm)');
    set(gca, 'view', [-50 20]);
    grid on;
end
figure(1); xl = xlim; figure(2); xlim(xl);
figure(1); yl = ylim; figure(2); ylim(yl);
figure(1); zl = zlim; figure(2); zlim(zl);

% convert to covariance, and display std & correlation
Sig = inv(V);
comp = sqrt(diag(Sig));
Corr = Sig./(comp*comp.');
disp(comp);
disp(Corr);

% draw the spatial maps
[LX, LY] = meshgrid(-.5:.5:7.5, -.5:.5:5.5);
for p = 1 : 3
    figure(p+2);
    GPR = fitrgp(location, dp(p, :)', 'FitMethod', 'none',  'KernelParameters', [2; sqrt(Sig(p, p)/2)]);
    LC = reshape(predict(GPR, [LX(:) LY(:)]), size(LX));
    surf(LX*5.2, -LY*5, LX*0, LC, 'EdgeColor', 'none', 'FaceColor', 'interp');
    view([0 90]);
    axis equal;
    axis off;
end


% fast discrete Gaussian calculation
function [mu, s] = intGaussian(m, lambda)
    global M LAM MU VA
    mu = zeros(size(m)); s = mu;
    range_lo = lambda < 1;
    range_hi = lambda > 80;
    range_mi = ~(range_lo | range_hi);
    res = m - round(m);
    mu(range_lo) = m(range_lo);
    s(range_lo) = 1./lambda(range_lo);
    mu(range_mi) = m(range_mi) + sign(res(range_mi)).*interp2(M, LAM, MU, abs(res(range_mi)), lambda(range_mi));
    s(range_mi) = interp2(M, LAM, VA, abs(res(range_mi)), lambda(range_mi))./lambda(range_mi);
    if any(range_hi)
	dmu = 1./(1 + exp(lambda(range_hi).*(0.5 - abs(res(range_hi)))));
	mu(range_hi) = m(range_hi) + sign(res(range_hi)).*dmu - res(range_hi);
	s(range_hi) = dmu.*(1-dmu);
    end
end
