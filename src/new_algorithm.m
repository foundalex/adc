function [sig_adc, delta_tilda] = new_algorithm(adc_input, s_to_subadc_int, s_after_subadc, sim_options);

    % save (sprintf(num2str(clock) + ".mat"));
    % load ('2026              3             18             15             41         46.311.mat');


    %% подсчет количества годных изделий
    Q = 2; % необходимый порядок ряда Тейлора
    wbeta = 0.88; % нормированная полоса пропускания
    tau = 0.01; % ско time skew
    M = 8; % количество каналов
    A = 1/(wbeta*pi) * (sqrt(2/3) * factorial(Q+1)/4096)^(1/(Q+1));
    P1_12 = A/(tau*(sqrt(2*(1-1/M))));
    erf_P1_12 = erf(P1_12);
    percent_erf = erf_P1_12*100;

    %% параметры моделирования
    N = 8192;
    hdc = 9;
    hdq1 = 25;
    hdq2 = 5;

    zfq1 = zeros(1,hdq1-1);
    zfq2 = zeros(1,hdq2-1);

    n = (0:1:sim_options.N-1);

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
    coeff = zeros(sim_options.M,1);

    del_proc0 = (hdc-1)/2;
    del_proc1 = (hdq1-1)/2;
    del_proc2 = (hdq2-1)/2;

    if sim_options.M == 4
        U = [2 -1 0 -1; -1 2 -1 0; 0 -1 2 -1; -1 0 -1 2];
    elseif sim_options.M == 2
        U = [2 -1; -1 2];
    end

    % находим псевдоинверсную матрицу
	pseudo_u = pinv(U);

    y1_mult = zeros(N*sim_options.M,1);

    %% 9 - ый порядок
    n_hdc_impz = -(hdc-1)/2:(hdc-1)/2;
    n_hdc = 0:hdc-1;

    w_blackman_hdc = 0.42 - 0.5 * cos(2*pi*n_hdc/(hdc-1)) + 0.08 * cos(4*pi*n_hdc/(hdc-1));
    for i = 1:hdc
        hdc_coeff(i) = ((-1)^n_hdc_impz(i))/n_hdc_impz(i); 
    end
    hdc_coeff((hdc+1)/2) = 0;

    Hdc = hdc_coeff .* w_blackman_hdc;
    Hd_hdc = dfilt.dfasymfir(Hdc);

    [x3_hdc, f3_hdc] = freqz(Hdc, 1);
    [x31_hdc, f31_hdc] = freqz(Hd_hdc.Numerator, 1);
    % figure(71); 
    % subplot(3,1,1)
    % plot(hdc_coeff)
    % title('Коэффициенты дифф.фильтра 9-го порядка')
    % xlabel('Номер коэффициента') 
    % ylabel('Значение коэффициента')
    % subplot(3,1,2)
    % plot(f3_hdc/pi, abs(x3_hdc), f3_hdc/pi, abs(x31_hdc));
    % title('АЧХ дифф.фильтра 9-го порядка')
    % xlabel('Нормированная частота (x Pi rad/sample)') 
    % ylabel('Амплитуда')
    % subplot(3,1,3)
    % plot(f3_hdc/pi, angle(x3_hdc), f3_hdc/pi, angle(x31_hdc));
    % title('ФЧХ дифф.фильтра 9-го порядка')
    % xlabel('Нормированная частота (x Pi rad/sample)') 
    % ylabel('Амплитуда')

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

    [x3, f3] = freqz(Hdq1, 1);
    [x31, f31] = freqz(Hdq11.Numerator);

    % figure(70); 
    % subplot(3,1,1)
    % plot(hdq1_coeff)
    % title('Коэффициенты дифф.фильтра 25-го порядка')
    % xlabel('Номер коэффициента') 
    % ylabel('Значение коэффициента')
    % subplot(3,1,2)
    % plot(f3/pi, abs(x3), f3/pi, abs(x31));
    % title('АЧХ дифф.фильтра 25-го порядка')
    % xlabel('Нормированная частота (x Pi rad/sample)') 
    % ylabel('Амплитуда')
    % subplot(3,1,3)
    % plot(f3/pi, angle(x3), f3/pi, angle(x31));
    % title('ФЧХ дифф.фильтра 25-го порядка')
    % xlabel('Нормированная частота (x Pi rad/sample)') 
    % ylabel('Амплитуда')

    %% 5 - ый порядок
    n_hdq2_impz = -(hdq2-1)/2:(hdq2-1)/2;
    n_hdq2 = 0:hdq2-1;
    w_blackman_hdq2 = 0.42 - 0.5 * cos(2*pi*n_hdq2/(hdq2-1)) + 0.08 * cos(4*pi*n_hdq2/(hdq2-1));

    for i = 1:hdq2
        hdq2_coeff(i) = ((-1)^n_hdq2_impz(i))/n_hdq2_impz(i); 
    end
    hdq2_coeff((hdq2+1)/2) = -1;

    % symmetric
    hdq2_coeff_symm = hdq2_coeff;
    hdq2_coeff_symm(4) = hdq2_coeff(2);
    hdq2_coeff_symm(5) = hdq2_coeff(1);

    % asymmetric
    
    

    Hdq2_symm = hdq2_coeff_symm .* w_blackman_hdq2;
    Hdq2_asymm = hdq2_coeff .* w_blackman_hdq2;
    % Hdq2 = dfilt.dfasymfir(Hdq2);
    % 
    % [x3_hdq2, f3_hdq2] = freqz(diff_filter_hdq2, 1,1024, 'whole', sim_options.Fs_sub_adc);
    [x31_hdq2, f31_hdq2] = freqz(Hdq2_symm, 1);
    [x31_hdq2_a, f31_hdq2_a] = freqz(Hdq2_asymm, 1);

    %%
    figure(70); 
    subplot(4,1,1)
    plot(Hdq2_symm)
    title('Коэффициенты дифф.фильтра 5-го порядка')
    xlabel('Номер коэффициента') 
    ylabel('Значение коэффициента')
    subplot(4,1,2)
    plot(Hdq2_asymm)
    title('Коэффициенты дифф.фильтра 5-го порядка')
    xlabel('Номер коэффициента') 
    ylabel('Значение коэффициента')
    subplot(4,1,3)
    plot(hdq1_coeff)
    title('Коэффициенты дифф.фильтра 25-го порядка')
    xlabel('Номер коэффициента') 
    ylabel('Значение коэффициента')
    subplot(4,1,4)
    plot(hdc_coeff)
    title('Коэффициенты дифф.фильтра 9-го порядка')
    xlabel('Номер коэффициента') 
    ylabel('Значение коэффициента')


    figure(71); 
    subplot(2,1,1)
    plot(f3_hdc/pi, abs(x3_hdc), f3/pi, abs(x3), f31_hdq2/pi, abs(x31_hdq2), f31_hdq2/pi, abs(x31_hdq2_a));
    title('АЧХ дифф.фильтров')
    xlabel('Нормированная частота (x Pi rad/sample)') 
    ylabel('Амплитуда')
    legend({'9-ый порядок', '25-ый порядок', '5-ый порядок'}, 'Location','northwest');
    subplot(2,1,2)
    plot(f31_hdq2/pi, angle(x3_hdc), f31_hdq2/pi, angle(x3), f31_hdq2/pi, angle(x31_hdq2), f31_hdq2/pi, angle(x31_hdq2_a));
    title('ФЧХ дифф.фильтров')
    xlabel('Нормированная частота (x Pi rad/sample)') 
    ylabel('Амплитуда')
    legend({'9-ый порядок', '25-ый порядок', '5-ый порядок'}, 'Location','northwest');


    %% Начало алгоритма
    % H(d,c)

    % y = filter(Hdc, double(s_after_subadc));
    y = filter(Hdc, 1, s_after_subadc);
    y_cut_s = y(del_proc0+1:end);

    % figure(55);
    % subplot(3,1,1)
    % plot([s_after_subadc(1:100), y_cut_s(1:100)])
    % subplot(3,1,2)
    % snr(s_after_subadc, sim_options.Fs_sub_adc*sim_options.M);
    % subplot(3,1,3)
    % snr(y_cut_s, sim_options.Fs_sub_adc*sim_options.M);
    %%
    
    % создаем массив каналов
    % ch_array = zeros(32768,1);
    % s_ind = 1;
    % e_ind = 4;
    % for r = 1:N
    %     ch_array(s_ind:e_ind) = (1:4);
    %     s_ind = s_ind + 4;
    %     e_ind = e_ind + 4;
    % end

    % buffer_hdq1 = zeros(1,hdq1);
    % buffer_hdq1_channel = zeros(1,hdq1);

    for j = 1:floor(length(y_cut_s)/(N*sim_options.M)) 

        % start_index = sim_options.M*N*(j-1)+1;
        % end_index = start_index+(N*sim_options.M)-1;

        start_index = sim_options.M*N*(j-1)+1;
        end_index = start_index+(N*sim_options.M)-1;

        for i = 1:sim_options.M
            y1_mult(i:sim_options.M:sim_options.M*N) = y_cut_s(start_index+i-1:sim_options.M:end_index) .* coeff(i,j);
        end
   
        % y1_mult(1:sim_options.M:sim_options.M*N) = y_cut_s(start_index:sim_options.M:end_index) .* coeff(4,j);
        % y1_mult(2:sim_options.M:sim_options.M*N) = y_cut_s(start_index+1:sim_options.M:end_index) .* coeff(1,j);
        % y1_mult(3:sim_options.M:sim_options.M*N) = y_cut_s(start_index+2:sim_options.M:end_index) .* coeff(2,j);
        % y1_mult(4:sim_options.M:sim_options.M*N) = y_cut_s(start_index+3:sim_options.M:end_index) .* coeff(3,j);

        y_tilda = s_after_subadc(start_index:end_index) - y1_mult;

        %% Hdq1
        
        % for n = 1:length(y_tilda)
        % 
        %     buffer_hdq1 = [(y_tilda(n)) buffer_hdq1(1:end-1)];
        %     buffer_hdq1_channel = [(ch_array(n)) buffer_hdq1_channel(1:end-1)];
        % 
        %     for i = 1:hdq1
        %         y_mult(i,n) =  Hdq1(i) * buffer_hdq1(i);
        %     end
        % 
        %     y_add(1,n) = y_mult(1,n) + y_mult(2,n);
        % 
        %     for i = 1:hdq1-2
        %         y_add(i+1,n) = y_add(i,n) + y_mult(i+2,n);
        %     end
        % end
        % 
        % y1_mf = (y_add(hdq1-1,:))';

        % y11  = filter(Hdq11, y_tilda);
        y_tilda = [y_tilda; zeros(del_proc1,1)];

        % figure(13)
        % plot(y_tilda1);


        [y1, zfq1] = filter(Hdq1, 1, y_tilda, zfq1);

        % if j == 1
            % y1_cut = [y1(del_proc1+1:end); zeros(del_proc1,1)];
            y1_cut = y1(del_proc1+1:end);
        % else
        %     y1_cut = y1;
        % end

        % figure(14)
        % plot(y1_cut);
    
        for i = 1:sim_options.M
            y1_cut(i:sim_options.M:sim_options.M*N) = y1_cut(i:sim_options.M:sim_options.M*N) .* coeff(i,j);
        end    

        
        % y1_cut(i:sim_options.M:sim_options.M*N) = y1_cut(i:sim_options.M:sim_options.M*N) .* coeff(4,j);
        % y1_cut(i:sim_options.M:sim_options.M*N) = y1_cut(i:sim_options.M:sim_options.M*N) .* coeff(1,j);
        % y1_cut(i:sim_options.M:sim_options.M*N) = y1_cut(i:sim_options.M:sim_options.M*N) .* coeff(2,j);
        % y1_cut(i:sim_options.M:sim_options.M*N) = y1_cut(i:sim_options.M:sim_options.M*N) .* coeff(3,j);

        %% Hdq2
        % [y2] = filter(Hdq2, y_tilda);

        y_tilda = [y_tilda; zeros(del_proc2,1)];

        [y2, zfq2] = filter(Hdq2_asymm, 1, y_tilda, zfq2);


        y2_cut = y2(del_proc2+1:end); % удаляем переходной процесс

        % figure(14)
        % plot(y2_cut);

        y2_cut1 = [y2_cut; zeros(del_proc2,1)];
        y2_cut = filter(Hdq2_symm, 1, y2_cut1);
        y2_cut = y2_cut(del_proc1+1:end-2);


        for i = 1:sim_options.M
            y2_cut(i:sim_options.M:sim_options.M*N) = y2_cut(i:sim_options.M:sim_options.M*N) .* ((coeff(i,j)^2)/2);
        end

        % y2_cut(i:sim_options.M:sim_options.M*N) = y2_cut(i:sim_options.M:sim_options.M*N) .* ((coeff(4,j)^2)/2);
        % y2_cut(i:sim_options.M:sim_options.M*N) = y2_cut(i:sim_options.M:sim_options.M*N) .* ((coeff(1,j)^2)/2);
        % y2_cut(i:sim_options.M:sim_options.M*N) = y2_cut(i:sim_options.M:sim_options.M*N) .* ((coeff(2,j)^2)/2);
        % y2_cut(i:sim_options.M:sim_options.M*N) = y2_cut(i:sim_options.M:sim_options.M*N) .* ((coeff(3,j)^2)/2);

        % скорректированный сигнал после дифф.фильтров
        x_tilda(:,j) = s_after_subadc(start_index:end_index) - y1_cut - y2_cut;

        %% производим корелляцию в окне и находим среднее значение
        for k = 1:sim_options.M-1
            cor(k,j) = sum(x_tilda(k:sim_options.M:sim_options.M*N,j).*x_tilda(k+1:sim_options.M:sim_options.M*N,j))/N;
        end
        cor(sim_options.M,j) = sum([0; x_tilda(sim_options.M:sim_options.M:sim_options.M*N-1,j)] .* x_tilda(1:sim_options.M:sim_options.M*N,j))/N;

        % находим ошибку среднего значения корреляций
        error_correlate(1,j) = cor(1,j) - cor(sim_options.M,j);
        for i = 2:sim_options.M
            error_correlate(i,j) = cor(i,j) - cor(i-1,j);
        end

		%% выход АЦП0 дифф. фильтра
        cor_M(j) = sum(x_tilda(2:sim_options.M:sim_options.M*N,j) .* y_cut_s(start_index:sim_options.M:end_index))/N;

		rate = nextpow2(cor_M(j));
		cor_M_floor = 2^rate;

		%%
		mult_u = pseudo_u * error_correlate(:,j);
		delta_tilda(:,j) = mult_u/cor_M_floor;
        coeff(:,j+1) = coeff(:,j) + 0.1 * delta_tilda(:,j);
        % coeff(:,j+1) = coeff(:,j);

        for i = 1:sim_options.M
            rms_delta(i,j) = delta_tilda(i,j); 
        end
    end

    % save (sprintf(num2str(clock) + ".mat"));
    % load ('2026              2             23             23             51         10.257.mat'); % f = 420MHz, ts = 0.05 

    start_index_s = 1;
    end_index_s = sim_options.M*N;

    sig_adc = zeros(sim_options.M*N*j,1);
    for i = 1:j
        sig_adc(start_index_s:end_index_s) = x_tilda(:,i);
        start_index_s = end_index_s + 1;
        end_index_s = start_index_s+(sim_options.M*N)-1;
    end

    figure(6)
    plot([s_to_subadc_int(end-100:end), sig_adc(end-100:end)]);

    %% Графики корреляций
    figure(7);
    for i = 1:sim_options.M
        subplot(sim_options.M,1,i)
        plot(rms_delta(i,:));
        title(['RMS delta' num2str(i)])
        xlabel('Номер значения') 
        ylabel('Значение') 
    end

    figure(8)
    subplot(sim_options.M+4,1,1)
    plot([s_to_subadc_int(end-250:end), s_after_subadc(end-250:end), sig_adc(end-250:end)]);
    title('Входной сигнал и выходной сигнал')
    xlabel('Номер отсчета') 
    ylabel('Значение отсчета') 
    legend({'Вход без ошибок','Вход с ошибками', 'Выход'},'Location','northeast');
    %% дельты
    for i = 1:sim_options.M
        subplot(sim_options.M+4,1,i+1)
        plot(delta_tilda(i,:));
        title(['Дельта суб-АЦП' num2str(i)])
        xlabel('Номер итерации') 
        ylabel('Значение дельты') 
    end

    subplot(sim_options.M+4,1,sim_options.M+2)
    snr(s_to_subadc_int*2^-11, sim_options.Fs_sub_adc*sim_options.M);
    subplot(sim_options.M+4,1,sim_options.M+3)
    snr(s_after_subadc*2^-11, sim_options.Fs_sub_adc*sim_options.M);
    subplot(sim_options.M+4,1,sim_options.M+4)
	snr(sig_adc*2^-11, sim_options.Fs_sub_adc*sim_options.M);

end