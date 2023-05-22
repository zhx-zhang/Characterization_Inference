% extract the peak locations and heights of the TM measurement from the selected list (xlist) of devices at the selected range of wavelength

load data/TM_data
range_s = lambda > 1540.241 & lambda < 1578.495;
l_s = lambda(range_s);
ng0 = 1.77;

xlist = [22:23 28:30 34:36 40:42 46:48];
TM_pk = cell(48, 1);

for n = xlist
    fprintf('starting: %d\n', n);
    R = radius(location(n, 1)+1)*1e3;
    b = data(range_s, n);
    base = smooth(lambda, data(:, n), 0.1, 'loess');
    b_dt = b - base(range_s);
    FSR = 1550^2/(ng0*2*pi*R);
    [pks, locs] = findpeaks(b_dt, l_s, 'MinPeakProminence', 10, 'MinPeakDistance', FSR*0.9);
    TM_pk{n} = struct('lambda', locs, 'val', pks);
end

save data/extract_TM TM_pk xlist
