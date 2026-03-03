function [sig_adc, delta_tilda] = new_algorithm(adc_input, s_to_subadc_int, s_after_subadc, sim_options);

    % load ('2026              3              3             17             38         30.605.mat');
    % save (sprintf(num2str(clock) + ".mat"));

    N = 8192;
    hdc = 9;
    hdq1 = 25;
    hdq2 = 5;

    zfq1 = zeros(1,hdq1-1);
    zfq2 = zeros(1,hdq2-1);

    n = (0:1:sim_options.N-1);
    % del_proc = ((sim_options.N-1)/2);

    % ws1 = 0.998;
    % ws2 = 0.001;
    % w_blackman = 0.42 - 0.5 * cos(2*pi*n/(sim_options.N-1)) + 0.08 * cos(4*pi*n/(sim_options.N-1));

    %%
    bb = firpm(hdc,[0 0.9],[0 0.9*pi],'d');
    Hd = dfilt.dfasymfir(bb);

    bb1 = firpm(hdq1,[0.1 0.9],[1 1],'differentiator');
    Hd1(1) = dfilt.dfasymfir(bb1);

    % figure(8);
    % subplot(1,1,1)
    % plot(Hd1(1,1).Numerator);

    Hd1s = dfilt.dfsymfir(bb1);

    bb2 = firpm(hdq2,[0.1 0.9],[1 1],'differentiator');
    Hd2 = dfilt.dfasymfir(bb2);
    % blo = firls(18,[0 0.45 0.55 1],[1 1 0 0],[100 1], 'differentiator');

    % [x, f] = freqz(Hd.Numerator, 1,1024, 'whole', sim_options.Fs_sub_adc);
    % [x1, f1] = freqz(Hd1(1).Numerator, 1,1024, 'whole', sim_options.Fs_sub_adc);
    % [xx1, ff1] = freqz(Hd1s.Numerator, 1,1024, 'whole', sim_options.Fs_sub_adc);
    % [x2, f2] = freqz(Hd2.Numerator, 1,1024, 'whole', sim_options.Fs_sub_adc);

    % fvtool(Hd,Hd1,Hd2);

    %%
    % figure(3);
    % subplot(3,1,1)
    % plot(f, abs(x));
    % subplot(3,1,2)
    % plot(f1, abs(x1));
    % subplot(3,1,3)
    % plot(f2, abs(x2));
    % % 
    % figure(4);
    % subplot(4,1,1)
    % plot(f, angle(x));
    % subplot(4,1,2)
    % plot(f1, angle(x1));
    % title('ФЧХ фильтра дифф.фильтра 25-го порядка assym')
    % xlabel('Частота') 
    % ylabel('Амплитуда') 
    % subplot(4,1,3)
    % plot(ff1, angle(xx1));
    % title('ФЧХ фильтра дифф.фильтра 25-го порядка symm')
    % xlabel('Частота') 
    % ylabel('Амплитуда')
    % subplot(4,1,4)
    % plot(f2, angle(x2));
    %%
    % figure(5);
    % subplot(2,1,1)
    % plot(f1, angle(x1));
    % title('ФЧХ фильтра дифф.фильтра 25-го порядка')
    % xlabel('Частота') 
    % ylabel('Амплитуда') 
    % subplot(2,1,2)
    % plot(Hd1.Numerator);
    % title('Импульсная характеристика')
    % xlabel('Номер отсчета') 
    % ylabel('Амплитуда') 

    %%
    coeff(1,1) = 0;
    coeff(2,1) = 0;
    coeff(3,1) = 0;
    coeff(4,1) = 0;

    del_proc0 = (hdc-1)/2;
    del_proc1 = (hdq1-1)/2;
    del_proc2 = (hdq2-1)/2;

    U = [2 -1 0 -1; -1 2 -1 0; 0 -1 2 -1; -1 0 -1 2];
    % находим псевдоинверсную матрицу
	pseudo_u = pinv(U);

    y1_mult = zeros(N*sim_options.M,1);

    %% 25 - ый порядок
    n1_impz = -(hdq1-1)/2:(hdq1-1)/2;
    n_hdq1 = 0:hdq1-1;

    w_blackman = 0.42 - 0.5 * cos(2*pi*n_hdq1/(hdq1-1)) + 0.08 * cos(4*pi*n_hdq1/(hdq1-1));
    for i = 1:hdq1
        hdq1_coeff(i) = ((-1)^n1_impz(i))/n1_impz(i); 
    end
    hdq1_coeff((hdq1+1)/2) = 0;

    % symmetric
    % hdq1_coeff(14) = hdq1_coeff(12);
    % hdq1_coeff(15) = hdq1_coeff(11);
    % hdq1_coeff(16) = hdq1_coeff(10);
    % hdq1_coeff(17) = hdq1_coeff(9);
    % hdq1_coeff(18) = hdq1_coeff(8);
    % hdq1_coeff(19) = hdq1_coeff(7);
    % hdq1_coeff(20) = hdq1_coeff(6);
    % hdq1_coeff(21) = hdq1_coeff(5);
    % hdq1_coeff(22) = hdq1_coeff(4);
    % hdq1_coeff(23) = hdq1_coeff(3);
    % hdq1_coeff(24) = hdq1_coeff(2);
    % hdq1_coeff(25) = hdq1_coeff(1);

    % a1 == a2 
    % for i = 1:N_taps
    %     a2(i) = cos(n11(i)*pi)/n11(i); 
    % end
    % a2(1) = 0;
    % a2(13) = 0;

    % figure(71); 
    % plot([a1,a2]);

    Hdq1 = hdq1_coeff .* w_blackman;

    Hdq11 = dfilt.dfasymfir(Hdq1);

    [x3, f3] = freqz(Hdq1, 1,1024, 'whole', 4*sim_options.Fs_sub_adc);
    [x31, f31] = freqz(Hdq11.Numerator, 1,1024, 'whole', 4*sim_options.Fs_sub_adc);

    figure(70); 
    subplot(3,1,1)
    plot(hdq1_coeff)
    subplot(3,1,2)
    plot(f3, abs(x3), f3, abs(x31));
    subplot(3,1,3)
    plot(f3, angle(x3), f3, angle(x31));

    %% 9 - ый порядок
    n_hdc_impz = -(hdc-1)/2:(hdc-1)/2;
    n_hdc = 0:hdc-1;

    w_blackman_hdc = 0.42 - 0.5 * cos(2*pi*n_hdc/(hdc-1)) + 0.08 * cos(4*pi*n_hdc/(hdc-1));
    for i = 1:hdc
        hdc_coeff(i) = ((-1)^n_hdc_impz(i))/n_hdc_impz(i); 
    end
    hdc_coeff((hdc+1)/2) = 0;

    Hdc = hdc_coeff .* w_blackman_hdc;
    % Hdc = dfilt.dfasymfir(Hdc);

    % 
    % [x3_hdc, f3_hdc] = freqz(diff_filter_hdc, 1,1024, 'whole', sim_options.Fs_sub_adc*4);
    % [x31_hdc, f31_hdc] = freqz(Hd_hdc.Numerator, 1,1024, 'whole', sim_options.Fs_sub_adc*4);
    % 
    % figure(71); 
    % subplot(3,1,1)
    % plot(hdc_coeff)
    % subplot(3,1,2)
    % plot(f3_hdc, abs(x3_hdc), f3_hdc, abs(x31_hdc), f3_hdc, abs(x));
    % subplot(3,1,3)
    % plot(f3_hdc, angle(x3_hdc), f3_hdc, angle(x31_hdc), f3_hdc, angle(x));
    % title('ФЧХ фильтра дифф.фильтра 9-го порядка')
    % xlabel('Частота') 
    % ylabel('Амплитуда')

    %% 5 - ый порядок
    n_hdq2_impz = -(hdq2-1)/2:(hdq2-1)/2;
    n_hdq2 = 0:hdq2-1;
    w_blackman_hdq2 = 0.42 - 0.5 * cos(2*pi*n_hdq2/(hdq2-1)) + 0.08 * cos(4*pi*n_hdq2/(hdq2-1));

    for i = 1:hdq2
        hdq2_coeff(i) = ((-1)^n_hdq2_impz(i))/n_hdq2_impz(i); 
    end
    hdq2_coeff((hdq2+1)/2) = 1;
    hdq2_coeff(4) = hdq2_coeff(2);
    hdq2_coeff(5) = hdq2_coeff(1);

    Hdq2 = hdq2_coeff .* w_blackman_hdq2;
    % Hdq2 = dfilt.dfasymfir(Hdq2);
    % 
    % [x3_hdq2, f3_hdq2] = freqz(diff_filter_hdq2, 1,1024, 'whole', sim_options.Fs_sub_adc);
    % [x31_hdq2, f31_hdq2] = freqz(Hd_hdq2.Numerator, 1,1024, 'whole', sim_options.Fs_sub_adc);
    % figure(72); 
    % subplot(3,1,1)
    % plot(hdq2_coeff)
    % subplot(3,1,2)
    % plot(f3_hdq2, abs(x3_hdq2), f3_hdq2, abs(x31_hdq2), f3_hdq2, abs(x2));
    % subplot(3,1,3)
    % plot(f3_hdc, angle(x3_hdq2), f3_hdq2, angle(x31_hdq2), f3_hdq2, angle(x2));
    % title('ФЧХ фильтра дифф.фильтра 5-го порядка')
    % xlabel('Частота') 
    % ylabel('Амплитуда')

    %%
    % figure(80); 
    % plot(f3_hdc, angle(x3_hdc), f3_hdc, angle(x3), f3_hdc, angle(x3_hdq2));

    %% VI samples intervention
    % in = 0;
    % row = 1;
    % start = 1;
    % 
    % while in < numel(s_after_subadc)-4
    %     for j = start:4
    %         in = in + 1;
    %         if in ~= 1 && mod(in-1,11) == 0
    %             j = 1;
    %             row = row + 1;
    %             sig(row,j) = s_after_subadc(in);
    %             start = 2;
    %         else
    %             sig(row,j) = s_after_subadc(in);
    %             start = 1;
    %         end
    %     end
    %     if mod(in-1,11) ~= 0
    %         row = row + 1;
    %     end
    % end
    % 
    % % удаление нулей из 4-го АЦП
    % odd_i = sig(1:3:end,4);
    % even_i = [sig(2:3:end,4)];
    % result = reshape([odd_i'; even_i'],1,[]);
    % 
    % adc_input_interv = zeros(numel(sig(:,1)), sim_options.M);
    % 
    % for i = 1:sim_options.M
    %     if i ~= sim_options.M
    %         adc_input_interv(:,i) = sig(:,i);
    %     else
    %         adc_input_interv(1:numel(result),i) = result;
    %     end
    % end

    %% H(d,c)

    % y = filter(Hdc, double(s_after_subadc));
    y = filter(Hdc, 1, s_after_subadc);
    y_cut_s = y(del_proc0+1:end);

    % figure(55);
    % subplot(2,1,1)
    % snr(double(adc_input(:,1)), sim_options.Fs_sub_adc);
    % subplot(2,1,2)
    % snr(y_cut(:,1), sim_options.Fs_sub_adc);
    %%
    % aa = zeros(32768,1);
    for j = 1:floor(length(y_cut_s)/(N*sim_options.M)) 

        start_index = sim_options.M*N*(j-1)+1;
        end_index = start_index+(N*sim_options.M)-1;
   
        y1_mult(1:sim_options.M:sim_options.M*N) = y_cut_s(start_index:sim_options.M:end_index) .* coeff(1,j);
        y1_mult(2:sim_options.M:sim_options.M*N) = y_cut_s(start_index+1:sim_options.M:end_index) .* coeff(2,j);
        y1_mult(3:sim_options.M:sim_options.M*N) = y_cut_s(start_index+2:sim_options.M:end_index) .* coeff(3,j);
        y1_mult(4:sim_options.M:sim_options.M*N) = y_cut_s(start_index+3:sim_options.M:end_index) .* coeff(4,j);

        y_tilda = s_after_subadc(start_index:end_index) - y1_mult;

        %% Hdq1
        % y11  = filter(Hdq11, y_tilda);

        [y1, zfq1] = filter(Hdq1, 1, y_tilda, zfq1);
        % if j == 1
            y1_cut = [y1(del_proc1+1:end); zeros(del_proc1,1)];
        % else
        %     y1_cut = y1;
        % end
        % 
        % figure(34);
        % plot([y1_cut(32669:32768), y1_cut(1:100)]);

        % aa = y1_cut;

        y1_cut(1:sim_options.M:sim_options.M*N) = y1_cut(1:sim_options.M:sim_options.M*N) .* coeff(1,j);
        y1_cut(2:sim_options.M:sim_options.M*N) = y1_cut(2:sim_options.M:sim_options.M*N) .* coeff(2,j);
        y1_cut(3:sim_options.M:sim_options.M*N) = y1_cut(3:sim_options.M:sim_options.M*N) .* coeff(3,j);
        y1_cut(4:sim_options.M:sim_options.M*N) = y1_cut(4:sim_options.M:sim_options.M*N) .* coeff(4,j);

        %% Hdq2
        % [y2] = filter(Hdq2, y_tilda);
        [y2, zfq2] = filter(Hdq2, 1, y_tilda, zfq2);
        % if j == 1
            y2_cut = [y2(del_proc2+1:end); zeros(del_proc2,1)];
        % else
        %     y2_cut = y2;
        % end

        y2_cut(1:sim_options.M:sim_options.M*N) = y2_cut(1:sim_options.M:sim_options.M*N) .* (coeff(1,j)^2)/2;
        y2_cut(2:sim_options.M:sim_options.M*N) = y2_cut(2:sim_options.M:sim_options.M*N) .* (coeff(2,j)^2)/2;
        y2_cut(3:sim_options.M:sim_options.M*N) = y2_cut(3:sim_options.M:sim_options.M*N) .* (coeff(3,j)^2)/2;
        y2_cut(4:sim_options.M:sim_options.M*N) = y2_cut(4:sim_options.M:sim_options.M*N) .* (coeff(4,j)^2)/2;

        x_tilda(:,j) = s_after_subadc(start_index:end_index) - y1_cut - y2_cut;

        % производим корелляцию в окне и находим среднее значение
        cor1(:,j) = sum(x_tilda(1:sim_options.M:sim_options.M*N,j).*x_tilda(2:sim_options.M:sim_options.M*N,j))/N;
        cor2(:,j) = sum(x_tilda(2:sim_options.M:sim_options.M*N,j).*x_tilda(3:sim_options.M:sim_options.M*N,j))/N;
        cor3(:,j) = sum(x_tilda(3:sim_options.M:sim_options.M*N,j).*x_tilda(4:sim_options.M:sim_options.M*N,j))/N;
        cor4(:,j) = sum([0; x_tilda(4:sim_options.M:sim_options.M*N-1,j)] .* x_tilda(1:sim_options.M:sim_options.M*N,j))/N;

        % находим ошибку среднего значения корреляций
		er(1,j) = cor1(:,j) - cor4(:,j);
		er(2,j) = cor2(:,j) - cor1(:,j);
		er(3,j) = cor3(:,j) - cor2(:,j);
		er(4,j) = cor4(:,j) - cor3(:,j);

		%% выход АЦП0 дифф. фильтра

        cor_M(j) = sum(x_tilda(2:sim_options.M:sim_options.M*N,j) .* y_cut_s(start_index:sim_options.M:end_index))/N;

		rate = nextpow2(cor_M(j));
		cor_M_floor = 2^rate;

		%%
		mult_u = pseudo_u * er(:,j);

		delta_tilda(:,j) = (mult_u/cor_M_floor);

        coeff(:,j+1) = coeff(:,j) + 0.2 * delta_tilda(:,j);

    end

    % save (sprintf(num2str(clock) + ".mat"));
    % load ('2026              2             23             23             51         10.257.mat'); % f = 420MHz, ts = 0.05 
 
    % выход АЦП
    % sig_adc = zeros(sim_options.M*length(x_tilda(:,1)),1);

    % for i = 1:numel(x_tilda(:,1))
    %     for j = start:4
    %         in = in + 1;
    %         if in ~= 1 && mod(in-1,11) == 0
    %             j = 1;
    %             row = row + 1;
    %             sig_adc(row) = x_tilda(i,j);
    %             start = 2;
    %         else
    %             row = row + 1;
    %             sig_adc(row) = x_tilda(i,j);
    %             start = 1;
    %         end
    %     end
    %     % if mod(in-1,11) ~= 0
    %     %     row = row + 1;
    %     % end
    % end

    % удаление нулей из 4-го АЦП

    % adc4 = x_tilda(1:363662,4);

    start_index_s = 1;
    end_index_s = sim_options.M*N;

    sig_adc = zeros(sim_options.M*N*j,1);
    for i = 1:j
        sig_adc(start_index_s:end_index_s) = x_tilda(:,i);
        start_index_s = end_index_s + 1;
        end_index_s = start_index_s+(sim_options.M*N)-1;
    end

    figure(6)
    subplot(8,1,1)
    plot([s_to_subadc_int(end-500:end), s_after_subadc(end-500:end), sig_adc(end-500:end)]);
    title('Входной сигнал и выходной сигнал')
    xlabel('Номер отсчета') 
    ylabel('Значение отсчета') 
    legend({'Вход','Выход'},'Location','northeast');
    %% дельты
    subplot(8,1,2)
    plot(delta_tilda(1,:));
    title('Дельта первого суб-АЦП')
    xlabel('Номер итерации') 
    ylabel('Значение дельты') 
    subplot(8,1,3)
    plot(delta_tilda(2,:));
    title('Дельта второго суб-АЦП')
    subplot(8,1,4)
    plot(delta_tilda(3,:));
    title('Дельта третьего суб-АЦП')
    subplot(8,1,5)
    plot(delta_tilda(4,:));
    title('Дельта четвертого суб-АЦП')
    subplot(8,1,6)
    snr(s_to_subadc_int, sim_options.Fs_sub_adc*4);
    subplot(8,1,7)
    snr(s_after_subadc, sim_options.Fs_sub_adc*4);
    subplot(8,1,8)
	snr(sig_adc, sim_options.Fs_sub_adc*4);

    %%
    figure(7);
    subplot(4,1,1)
    plot(er(1,:));
    subplot(4,1,2)
    plot(er(2,:));
    subplot(4,1,3)
    plot(er(3,:));
    subplot(4,1,4)
    plot(er(4,:));

end