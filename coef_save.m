% calculate neff and ng from extracted peak features
% combine with coupling & loss from coupling_loss_assign.m as the final extraction data
% bend.mat: simulation data from MODE that contains sensitivity of neff and ng.
% (Note: radius 20 um not included in the simulation data)

load data/bend
lc = lambda(6)*1e9;
NEFF0 = squeeze(real(neff(6, :, 1))) - 0.006;

load data/extract_raw
load data/extract_coupling_loss

N_S = length(L);

neff = cell(N_S, 1);
ng = cell(N_S, 1);

for n = 1 : N_S
    extract = extract_info{n};
    Err = error_extract{n};
    R = radius(location(n, 1) + 1)*1e3;

    if isnan(extract)
	continue;
    end

    % find the relative mode number from fitted ng/FSR
    neff0 = NEFF0(location(n, 1));
    N0 = neff0*2*pi*R/lc;
    f = extract(:, 1); err_f = Err(:, 1);
    FSR = f(2:end) - f(1:end-1);
    err_FSR = sqrt(err_f(2:end).^2 + err_f(1:end-1).^2);
    Fc = (f(1:end-1) + f(2:end))/2;
    err_fc = err_FSR/2;
    NG = Fc.^2./FSR/(2*pi*R);
    err_ng = NG.*sqrt(4*(err_fc./Fc).^2 + (err_FSR./FSR).^2);
    xrange = abs(NG - 1.77) > .5;
    w1 = 1./ err_ng;
    w1(xrange) = 0;
    Xc = Fc.^(0:1);
    b = (w1.*Xc)\(w1.*NG);
    ng_fit = f.^(0:1)*b;
    FSR_fit = f.^2./ng_fit/(2*pi*R);
    f_pred = f(end);
    while f_pred(end) > f(1)
	FSR_pred = interp1(f, FSR_fit, f_pred(end));
	f_next = f_pred(end) - FSR_pred;
	f_pred = [f_pred; f_next];
    end
    [~, J] = min(abs(f' - f_pred));

    % use the value at 1574 nm to calibrate the absolute mode number
    f_coded = (f-1550)/50;
    X = f_coded.^(0:1);
    w = 1./ err_f;
    b = (w.*X)\(w.*J');
    r_fit = X*b - J';
    S = (w.*r_fit)'*(w.*r_fit);
    M = S/(length(f) - 2)*inv(X'*diag(w.^2)*X);
    lc_coded = (lc-1550)/50;
    Xc = lc_coded.^(0:1);
    fr = Xc*b;
    e_fr = sqrt(Xc*M*Xc.');
    dN = N0 - fr;
    N = fr + round(dN);
    v_neff = lc*N/(2*pi*R);
    d_neff = lc/(2*pi*R);
    e_neff = e_fr*lc/(2*pi*R);
    neff{n} = struct('lambda', lc, 'val', v_neff, 'gap', d_neff, 'err', e_neff);
    ng{n} = struct('lambda', Fc(~xrange), 'val', NG(~xrange), 'err', err_ng(~xrange));
end

R = 20*(location(:, 1)+1);
g = 0.6 + location(:, 2)*0.2;
xrange = cellfun(@isempty, neff);
neff(xrange) = []; ng(xrange) = []; L(xrange) = []; loga(xrange) = [];
R(xrange) = []; g(xrange) = []; location(xrange, :) = [];
save data/extract_final neff ng L loga location R g
