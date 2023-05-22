% extract the resonance information (peak_info) from fitting of the measurement, 
% and partially transform to physical properties (extract_info). 
% estimated errors are also calculated (error_peak), transformed (error_extract).
% baseline.mat: a pre-computed baseline response (y_sm) from RRs with no observed resonances
% peak_info / error_peak columns: location, height, width
% extract_info / error_extract columns: location, ng (fitted), coupling, loss 
% (the last two can be exchanged)

load data/ring_data
load data/baseline
range_s = lambda > 1540.241 & lambda < 1609.545;
period = 0.25;
l_s = lambda(range_s);
N_S = size(data, 2);
peak_info = cell(N_S, 1);
error_peak = cell(N_S, 1);
extract_info = cell(N_S, 1);
error_extract = cell(N_S, 1);
F_lorentz = @(x, lambda) x(1) - x(2)./(1 + ((lambda - x(3))/x(4)).^2);
opts = optimset('Display','off');
ng0 = 1.77;

for n = 1 : N_S
    R = radius(location(n, 1)+1)*1e3;
    y = data(range_s, n);
    y_dt = y - y_sm(range_s);
    y_dt = 10.^(y_dt/10);
    FSR = 1550^2/(ng0*2*pi*R);
    if R > 40e3
	tol = .08;
    else
	tol = .1;
    end
    % find resonances
    [~, locs, w, p] = findpeaks(-y_dt, l_s, 'MinPeakProminence', tol, 'MinPeakDistance', FSR*0.9);
    if length(locs) < 3
	peak_info{n} = NaN;
	extract_info{n} = NaN;
	error_peak{n} = NaN;
	error_extract{n} = NaN;
	continue;
    end
    N_P = length(locs);
    info = zeros(N_P, 3);
    d_info = zeros(N_P, 3);
    % feature extraction
    for peak = 1 : N_P
	ind = N_P + 1 - peak;
	FR_r(n, peak) = locs(ind);
	w_r(n, peak) = w(ind);
	E_r(n, peak) = p(ind);
	% fitting
	w_bound = 3*w(peak) + period;
	range_fit = abs(l_s - locs(peak)) <= w_bound;
	l_fit = l_s(range_fit);
	y_fit = y_dt(range_fit);
	max_l = l_fit(end); min_l = l_fit(1);
	scale = (max_l - min_l)/2;
	l_center = (max_l + min_l)/2;
	l_coded = (l_fit - l_center)/scale;
	x0 = [0 p(ind) 0 w(ind)/2/scale];
	[x, resnorm, ~, ~, ~, ~, J] = lsqcurvefit(F_lorentz, x0, l_coded, y_fit, [], [], opts);
	dx = sqrt(resnorm / (length(l_coded) - 4) * diag(inv(J'*J)));
	info(peak, 1) = l_center + x(3)*scale;
	info(peak, 2) = x(4)*scale*2;
	d_info(peak, 1) = dx(3)*scale;
	d_info(peak, 2) = dx(4)*scale*2;

	y_fit = 10*log10(y_fit);
	x0 = [0 max(-y_fit) x(3) x(4)];
	[x, resnorm, ~, ~, ~, ~, J] = lsqcurvefit(F_lorentz, x0, l_coded, y_fit, [], [], opts);
	dx = sqrt(resnorm / (length(l_coded) - 4) * diag(inv(J'*J)));
	info(peak, 3) = x(2);
	d_info(peak, 3) = dx(2);
    end
    pos = info(:, 3) > 1e-4 & info(:, 2) > 1e-4;
    info = info(pos, :);
    d_info = d_info(pos, :);
    peak_info{n} = info;
    error_peak{n} = d_info;

    % property transform
    FSR = info(2:end, 1) - info(1:end-1, 1);
    err_FSR = sqrt(d_info(2:end, 1).^2 + d_info(1:end-1, 1).^2);
    Fc = (info(1:end-1, 1) + info(2:end, 1))/2;
    err_fc = err_FSR/2;
    ng = Fc.^2./FSR/(2*pi*R);
    err_ng = ng.*sqrt(4*(err_fc./Fc).^2 + (err_FSR./FSR).^2);
    % figure(1); plot(Fc, ng);
    % fitting
    % w = (1 - 10.^(-info(:, 3)/10)).^2./info(:, 2);
    xrange = abs(ng - 1.77) > .5;
    % w1 = (w(1:end-1)+w(2:end))/2;
    w1 = 1./ err_ng;
    w1(xrange) = 0;
    Xc = Fc.^(0:1);
    X = info(:, 1).^(0:1);
    b = (w1.*Xc)\(w1.*ng);
    ng_fit = X*b;
    r_fit = Xc*b - ng;
    S = (w1.*r_fit)'*(w1.*r_fit);
    M = S/(length(ng) - 2) * inv(Xc'*diag(w1.^2)*Xc);
    err_ng_fit = sqrt(diag(X*M*X'));
    
    FSR_fit = info(:, 1).^2./ng_fit/(2*pi*R);
    err_FSR_fit = FSR_fit.*sqrt(4*(d_info(:, 1)./info(:, 1)).^2 + (err_ng_fit./ng_fit).^2);
    FN = FSR_fit./info(:, 2);
    pFN = pi./FN;
    err_F = pFN.* sqrt((err_FSR_fit./FSR_fit).^2 + (d_info(:, 2)./info(:, 2)).^2);
    ER = 10.^(info(:, 3)/10);
    iER = max(1./ER - 1/400, 0);
    err_E = log(10)/20 * iER.* d_info(:, 3);
    A = cos(pFN)./(1+sin(pFN));
    B = 1 - (1-cos(pFN))./(1+cos(pFN)).*iER;
    alpha = sqrt(A./B) + sqrt(A./B - A);
    t = sqrt(A./B) - sqrt(A./B - A);

    dir_F = (1 + sqrt(iER)).*err_F;
    dir_E = pFN.*err_E;
    err_t = sqrt(dir_F.^2 + dir_E.^2 + 1.3*dir_F.*dir_E)/2;
    dir_F2 = (1 - sqrt(iER)).*err_F;
    err_a = sqrt(dir_F2.^2 + dir_E.^2 - 1.3*dir_F2.*dir_E)/2;
    extract = [info(:, 1) ng_fit alpha t];
    extract_info{n} = extract;
    error_extract{n} = [d_info(:, 1) err_ng_fit err_a err_t];
    % figure(2); plot(info(:, 1), [alpha t]);
    % g = gap(location(n, 2)+1);
    % title(['R = '  num2str(R/1e3)  ', g = '  num2str(g)]);
    % xlabel('wavelength (nm)');
    % grid on;
    % pause;
end

save data/extract_raw peak_info error_peak extract_info error_extract location radius gap
