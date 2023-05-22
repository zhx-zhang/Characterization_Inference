% nominal fitting design and sensitivity model based on simulation data: bend.mat (waveguide), coupling.mat, coupling_base.mat (coupler)
% also scale the error estimation of the extracted values based on wavelength response

load data/bend
% neff: nominal fitting design and sensitivity
lc = lambda(6);
neff = real(squeeze(neff(6, :, :)));
a_neff = neff(:, 1);
a_neff = a_neff - mean(a_neff); a_neff = a_neff/sqrt(var(a_neff));
s_neff = (neff(:, 2:2:6) - neff(:, 3:2:7))/2./[65 10 3.5];

% ng: nominal fitting design and sensitivity
dn = (ng(:, :, 2:2:6) - ng(:, :, 3:2:7))/2./permute([65 10 3.5], [1 3 2]);
l_max = max(lambda); l_min = min(lambda);
scale_n = l_max - l_min;
l_n_c = (l_max + l_min)/2;
l_coded = (lambda - l_n_c)/scale_n*2;
lc_coded = l_coded(6);
X_lam = l_coded.^(0:1);
coef = X_lam \ reshape(dn, size(dn, 1), []);
c_dn = reshape(coef, size(coef, 1), size(dn, 2), []);
c_n = X_lam \ ng(:, :, 1);
c_n = c_n - mean(c_n, 2); c_n = c_n./sqrt(var(c_n, [], 2));

% L: sensitivity
load data/coupling
load data/coupling_base % the simulation of a single waveguide but with other setting exactly the same as in coupler simulation
Tb = permute(Tb, [1 3 2]); % use the base value to eliminate numerical error
t = T./Tb;
L = log(acos(abs(t)));
L_norm = mean(L, 3);
dL = (L(:, :, 1:2:5) - L(:, :, 2:2:6))/2./permute([65 10 3.5], [1 3 2]);
l_max = max(lambda); l_min = min(lambda);
scale_t = l_max - l_min;
l_t_c = (l_max + l_min)/2;
l_coded = (lambda - l_t_c)/scale_t*2;
X_lam = l_coded.^(0:2);
coef = X_lam \ reshape(dL, size(dL, 1), []);
coef = reshape(coef, size(coef, 1), size(dL, 2), []);
XR = x2fx(DoE./[100 1], 'quadratic');
b_dt = XR \ reshape(permute(coef, [2 1 3]), size(dL, 2), []);

% L: nominal fitting design
load data/extract_final
A = cell(size(R)); S = cell(size(R)); y = cell(size(R));
e_y = cell(size(R)); d_y = cell(size(R));
X_full = x2fx([R/100 g], 'quadratic');
c_dt = permute(reshape(X_full*b_dt, size(X_full, 1), size(coef, 1), []), [2 1 3]);

% error estimation adjustment & combining matrices
for n = 1 : length(R)
    % neff
    R_index = location(n, 1);
    se = s_neff(R_index, :);
    ae = a_neff(R_index).^(0:1);

    % ng
    l_n_coded = (ng{n}.lambda*1e-9 - l_n_c)/scale_n*2;
    X_lam_n = l_n_coded.^(0:1);
    c_sn = squeeze(c_dn(:, R_index, :));
    sn = X_lam_n*c_sn;

    c_an = c_n(:, R_index);
    an = X_lam_n * [eye(length(c_an)) diag(c_an)];

    w = 1./ng{n}.err;
    b = (w.*X_lam_n) \ (w.*ng{n}.val);
    const_n = norm(w.*(X_lam_n*b - ng{n}.val))/sqrt(size(X_lam_n)*[1; -1]);
    ng{n}.err = const_n*ng{n}.err;

    % L
    l_t_coded = (L{n}.lambda*1e-9 - l_t_c)/scale_t*2;
    X_lam_t = l_t_coded.^(0:2);
    c_st = squeeze(c_dt(:, n, :));
    st = X_lam_t*c_st;

    c_at = [X_full(n, :) (R(n)/100).^(3:4)];
    at = X_lam_t * kron(c_at, eye(size(X_lam_t, 2)));
    
    w = 1./L{n}.err;
    b = (w.*X_lam_t) \ (w.*L{n}.val);
    const_t = norm(w.*(X_lam_t*b - L{n}.val))/sqrt(size(X_lam_t)*[1; -1]);
    L{n}.err = const_t*L{n}.err;

    S{n} = [se, se*0; sn, sn*0; st, st*R(n)];
    A{n} = blkdiag(ae, an, at);
    y{n} = [neff{n}.val; ng{n}.val; L{n}.val];
    e_y{n} = [neff{n}.err; ng{n}.err; L{n}.err];
    d_y{n} = [neff{n}.gap; l_n_coded*0; l_t_coded*0];
end

save data/inference_input S A y e_y d_y location R g
