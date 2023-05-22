% assign correct coupling and loss coefficients based on the wavelength response and comparing across devices. details/derivation of the method are included in my PhD thesis: 
% Adjoint Methods and Inverse Modeling for Process Variation Analysis in Silicon Photonics,
% which should be published by my school soon.
% output will also transform the coupling coefficient into the L-coefficient discussed in the paper

load data/extract_raw
load data/extract_TM
save_gif = false;

L = cell(48, 1);
loga = cell(48, 1);

for R_start = 0:6:42

    % reading data and normalization
    Val = cell(6, 1);
    Err = cell(6, 1);
    Prob = cell(6, 1);
    x = cell(6, 1);

    R = radius(R_start/6+1);
    filename = ['figure/R', num2str(R), '.gif'];
    h = figure(1);
    for n = 1:6
        extract = extract_info{n+R_start};
        err = error_extract{n+R_start};
        if ismember(n+R_start, xlist)
            xindex = TM_pk{n+R_start}.val > 10;
    	    [~, xrange] = min(abs(TM_pk{n+R_start}.lambda(xindex)' - extract(:, 1)));
       	    extract(xrange, :) = [];
    	    err(xrange, :) = [];
        end
        if ~isnan(extract)
       	    Val{n} = extract(:, 3:4);
    	    Err{n} = err(:, 3:4);
       	    Prob{n} = [0.55 0.45] + Val{n}*0;
    	    x{n} = (extract(:, 1)-1550)/50;
        end
    end

    % model: find common loss
    x_exp = cat(1, x{:});
    if length(x_exp) == 0
    	continue;
    end
    X_exp = x_exp.^(0:2);

    for step = 1 : 20
        Val_exp = cat(1, Val{:});
        Prob_exp = cat(1, Prob{:});
        Val_avg = sum(Val_exp.*(1-Prob_exp), 2);
    	Prob_avg = mean(Prob_exp);
        Val_var = Prob_exp(:, 1).*Prob_exp(:, 2).*(Val_exp(:, 1) - Val_exp(:, 2)).^2;
        b_a = X_exp \ Val_avg;
        r_fit = X_exp*b_a - Val_avg;
        e_a = r_fit'*r_fit/(length(Val_avg) - 0) + mean(Val_var);

        for n = 1:6
            if length(x{n}) < 1
                hs = subplot(2, 3, n);
                delete(hs);
                continue;
            end
            X = x{n}.^(0:2);
    	    X_t = x{n}.^(0:1);

    	    if step > 3
		% extra model: linear coupling
		Val_avg = sum(Val{n}.* Prob{n}, 2);
                Val_var = Prob{n}(:, 1).*Prob{n}(:, 2).*(Val{n}(:, 1) - Val{n}(:, 2)).^2;
		b_t = X_t \ Val_avg;
		fit_t = X_t*b_t;
		r_fit = fit_t - Val_avg;

		e_t = r_fit'*r_fit/(length(Val_avg) - 0) + mean(Val_var);
    	    end

            fit_a = X*b_a;
            dLP = ((Val{n}(:, 2) - fit_a).^2 - (Val{n}(:, 1) - fit_a).^2)./e_a;
    	    if step > 3
		dLP2 = ((Val{n}(:, 1) - fit_t).^2 - (Val{n}(:, 2) - fit_t).^2)./e_t;
		dLP = dLP + dLP2;
    	    end
    	    p_norm = Prob_avg(1)./ (Prob_avg(2)*exp(dLP) + Prob_avg(1));
            Prob{n}(:, 1) = p_norm;
            Prob{n}(:, 2) = 1 - p_norm;

            ind_sw = p_norm < 0.5;
            Val{n}(ind_sw, 1:2) = Val{n}(ind_sw, [2 1]);
            Err{n}(ind_sw, 1:2) = Err{n}(ind_sw, [2 1]);
            Prob{n}(ind_sw, 1:2) = Prob{n}(ind_sw, [2 1]);
            subplot(2, 3, n);
            % errorbar(x{n} * 50 + 1550, Val{n}(:, 1), Err{n}(:, 1));
    	    plot(x{n} * 50 + 1550, Val{n}, '-o');
            hold on;
            % errorbar(x{n} * 50 + 1550, Val{n}(:, 2), Err{n}(:, 2));
            errorbar(x{n} * 50 + 1550, fit_a, sqrt(e_a)+x{n}*0, '--');
    	    if step > 3
		errorbar(x{n} * 50 + 1550, fit_t, sqrt(e_t)+x{n}*0, '--');
    	    end
    	    hold off;
            xlim([min(x{n}) max(x{n})] * 50 + 1550);
    	    % switch the assignment just for the plot temporarily
            Val{n}(ind_sw, 1:2) = Val{n}(ind_sw, [2 1]);
            Err{n}(ind_sw, 1:2) = Err{n}(ind_sw, [2 1]);
            Prob{n}(ind_sw, 1:2) = Prob{n}(ind_sw, [2 1]);
        end

        sgtitle(['R = ', num2str(R), ' um']);
        if save_gif
            drawnow;
            frame = getframe(h);
            im = frame2im(frame);
            [imind, cm] = rgb2ind(im, 256);
            if step == 1
                imwrite(imind, cm, filename, 'gif', 'Loopcount', inf);
            else
                imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append');
            end
    	else
    	    pause(0.1);
        end
    end

    for n = 1 : 6
    	if length(x{n}) < 1
    	    continue;
    	end
    	ind_sw = Prob{n}(:, 1) < 0.5;
       	Val{n}(ind_sw, 1:2) = Val{n}(ind_sw, [2 1]);
    	Err{n}(ind_sw, 1:2) = Err{n}(ind_sw, [2 1]);
       	Prob{n}(ind_sw, 1:2) = Prob{n}(ind_sw, [2 1]);
    	Fc = x{n} * 50 + 1550;
    	acost = acos(Val{n}(:, 1)); v_L = log(acost);
    	e_L = Err{n}(:, 1)./sin(acost)./acost;
    	L{R_start+n} = struct('lambda', Fc, 'val', v_L, 'err', e_L);

    	Ac = Val{n}(:, 2);
    	v_loga = -20*log10(Ac)/(2*pi*R*1e-4);
    	e_loga = 20*Err{n}(:, 2)./Ac/(log(10)*2*pi*R*1e-4);
    	loga{R_start+n} = struct('lambda', Fc, 'val', v_loga, 'err', e_loga);
    end
    if save_gif
	save(['data_R', num2str(R), '.mat'], 'Val', 'Err', 'Prob', 'x');
    end
end

save data/extract_coupling_loss L loga
