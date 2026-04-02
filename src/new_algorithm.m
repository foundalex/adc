function [sig_adc, delta_tilda] = new_algorithm(adc_input, s_to_subadc_int, s_after_subadc, sim_options);
    %% подсчет количества годных изделий
    Q = 2; % необходимый порядок ряда Тейлора
    wbeta = 0.88; % нормированная полоса пропускания
    tau = 0.01; % ско time skew
    M = 8; % количество каналов АЦП
    A = 1/(wbeta*pi) * (sqrt(2/3) * factorial(Q+1)/4096)^(1/(Q+1));
    P1_12 = A/(tau*(sqrt(2*(1-1/M))));
    erf_P1_12 = erf(P1_12);
    percent_erf = erf_P1_12*100;

    %% параметры моделирования
    N = 8192; % количество отсчетов для каждого канала
    hdc = 9; % порядок фильтра Hdc
    hdq1 = 25; % порядок фильтра Hdq1
    hdq2 = 5; % порядок фильтра Hdq2
    mu = 0.1; % шаг LMS

    del_proc0 = (hdc-1)/2; % переходной процесс для фильтра Hdc
    del_proc1 = (hdq1-1)/2; % переходной процесс для фильтра Hdq1
    del_proc2 = (hdq2-1)/2; % переходной процесс для фильтра Hdq2

    % матрица U для 4-ех и для 2-ух каналов
    if sim_options.M == 4
        U = [2 -1 0 -1; -1 2 -1 0; 0 -1 2 -1; -1 0 -1 2];
    elseif sim_options.M == 2
        U = [2 -1; -1 2];
    end

    % находим псевдоинверсную матрицу от U
	pseudo_u = pinv(U);

    coeff = zeros(sim_options.M,1); % начальные значения дельт
    y1_mult = zeros(N*sim_options.M,1); 

    %% Синтез фильтра 9-го порядка
    n_hdc_impz = -(hdc-1)/2:(hdc-1)/2;
    n_hdc = 0:hdc-1; % кол-во отсчетов для окна Блэкмена

    w_blackman_hdc = 0.42 - 0.5 * cos(2*pi*n_hdc/(hdc-1)) + 0.08 * cos(4*pi*n_hdc/(hdc-1)); % окно Блэкмена
    for i = 1:hdc
        hdc_coeff(i) = ((-1)^n_hdc_impz(i))/n_hdc_impz(i); % ИХ дифференциального фильтра
    end
    hdc_coeff((hdc+1)/2) = 0;

    Hdc = hdc_coeff .* w_blackman_hdc; % умножение ИХ на окно Блэкмена

    [x3_hdc, f3_hdc] = freqz(Hdc, 1); % ЧХ фильтра

    %% Синтез фильтра 25-го порядка
    n1_impz = -(hdq1-1)/2:(hdq1-1)/2;
    n_hdq1 = 0:hdq1-1;

    w_blackman = 0.42 - 0.5 * cos(2*pi*n_hdq1/(hdq1-1)) + 0.08 * cos(4*pi*n_hdq1/(hdq1-1));
    for i = 1:hdq1
        hdq1_coeff(i) = ((-1)^n1_impz(i))/n1_impz(i); 
    end
    hdq1_coeff((hdq1+1)/2) = 0;

    Hdq1 = hdq1_coeff .* w_blackman;

    [x3, f3] = freqz(Hdq1, 1); % ЧХ фильтра

    %% Синтез фильтра 5-го порядка
    n_hdq2_impz = -(hdq2-1)/2:(hdq2-1)/2;
    n_hdq2 = 0:hdq2-1;
    w_blackman_hdq2 = 0.42 - 0.5 * cos(2*pi*n_hdq2/(hdq2-1)) + 0.08 * cos(4*pi*n_hdq2/(hdq2-1));

    for i = 1:hdq2
        hdq2_coeff(i) = ((-1)^n_hdq2_impz(i))/n_hdq2_impz(i); 
    end
    hdq2_coeff((hdq2+1)/2) = -1;

    % делаем ИХ симметричной [1] стр. 8
    hdq2_coeff_symm = hdq2_coeff;
    hdq2_coeff_symm(4) = hdq2_coeff(2);
    hdq2_coeff_symm(5) = hdq2_coeff(1);
    
    Hdq2_symm = hdq2_coeff_symm .* w_blackman_hdq2; % симметричная ИХ
    Hdq2_asymm = hdq2_coeff .* w_blackman_hdq2; % асимметричная ИХ

    [x31_hdq2, f31_hdq2] = freqz(Hdq2_symm, 1); % ЧХ симметричной ИХ
    [x31_hdq2_a, f31_hdq2_a] = freqz(Hdq2_asymm, 1); % ЧХ асимметричной ИХ

    %% Графики ИХ фильтров
    figure(70); 
    subplot(3,1,1)
    plot(Hdq2_symm)
    title('Коэффициенты дифф.фильтра 5-го порядка')
    xlabel('Номер коэффициента') 
    ylabel('Значение коэффициента')
    subplot(3,1,2)
    plot(hdq1_coeff)
    title('Коэффициенты дифф.фильтра 25-го порядка')
    xlabel('Номер коэффициента') 
    ylabel('Значение коэффициента')
    subplot(3,1,3)
    plot(hdc_coeff)
    title('Коэффициенты дифф.фильтра 9-го порядка')
    xlabel('Номер коэффициента') 
    ylabel('Значение коэффициента')

    %% Графики ЧХ фильтров
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


    %% ******************************************************************** Начало алгоритма калибровки

    % фильтр Hd,c - для coarse correction
    y = filter(Hdc, 1, s_after_subadc);
    % удаление переходного процесса фильтра
    y_cut_s = y(del_proc0+1:end);
    
    % работа с каждыми N*M отсчетами
    for j = 1:floor(length(y_cut_s)/(N*sim_options.M)) 

        % начало каждого N*M пакета 
        start_index = sim_options.M*N*(j-1)+1;
        % конец каждого N*M пакета 
        end_index = start_index+(N*sim_options.M)-1;

        % умножение выходных отсчетов на соответствующую рассчитаную дельту
        % для каждого канала
        for i = 1:sim_options.M
            y1_mult(i:sim_options.M:sim_options.M*N) = y_cut_s(start_index+i-1:sim_options.M:end_index) .* coeff(i,j);
        end

        % Coarse correction
        y_tilda = s_after_subadc(start_index:end_index) - y1_mult;

        %% Fine correction

        % добавление нулей для выравнивания количества отсчетов на выходе
        % параллельных фильтров
        y_tilda = [y_tilda; zeros(del_proc1,1)];

        % фильтр Hd q1
        y1 = filter(Hdq1, 1, y_tilda);
        % удаление переходного процесса фильтра
        y1_cut = y1(del_proc1+1:end);

        % умножение выходных отсчетов на соответствующую рассчитаную дельту
        % для каждого канала
        for i = 1:sim_options.M
            y1_cut(i:sim_options.M:sim_options.M*N) = y1_cut(i:sim_options.M:sim_options.M*N) .* coeff(i,j);
        end    

        % фильтр Hd q2_1 первая производная
        y2 = filter(Hdq2_asymm, 1, y_tilda);
        % удаление переходного процесса фильтра
        y2_cut = y2(del_proc2+1:end); 

        % фильтр Hd q2_2 вторая производная
        y2_cut = filter(Hdq2_symm, 1, y2_cut);
        % удаление переходного процесса фильтра
        y2_cut = y2_cut(del_proc2+1:end-8);

        % умножение выходных отсчетов на соответствующую рассчитаную дельту
        % для каждого канала
        for i = 1:sim_options.M
            y2_cut(i:sim_options.M:sim_options.M*N) = y2_cut(i:sim_options.M:sim_options.M*N) .* ((coeff(i,j)^2)/2);
        end

        % скорректированный сигнал после дифф.фильтров
        x_tilda(:,j) = s_after_subadc(start_index:end_index) - y1_cut - y2_cut;

        %% производим корелляцию в окне между каналами и находим среднее значение
        for k = 1:sim_options.M-1
            cor(k,j) = sum(x_tilda(k:sim_options.M:end-1,j).*x_tilda(k+1:sim_options.M:end,j))/N;
        end
        % последний канал задерживаем на 1 отсчет
        cor(sim_options.M,j) = sum([0; x_tilda(sim_options.M:sim_options.M:end-1,j)] .* x_tilda(1:sim_options.M:end,j))/N;

        % находим ошибку среднего значения корреляций
        error_correlate(1,j) = cor(1,j) - cor(sim_options.M,j);
        for i = 2:sim_options.M
            error_correlate(i,j) = cor(i,j) - cor(i-1,j);
        end

		% поиск производной [1] (15)
        cor_M(j) = sum(x_tilda(2:sim_options.M:sim_options.M*N,j) .* y_cut_s(start_index:sim_options.M:end_index))/N;
		rate = nextpow2(cor_M(j));
		cor_M_floor = 2^rate;

		% поиск вектора дельты [1] (16)
		mult_u = pseudo_u * error_correlate(:,j);
		delta_tilda(:,j) = mult_u/cor_M_floor;

        % LMS алгоритм
        coeff(:,j+1) = coeff(:,j) + mu * delta_tilda(:,j);

        % рассчитываем RMS дельты для графика
        for i = 1:sim_options.M
            rms_delta(i,j) = delta_tilda(i,j); 
        end
    end

    % собираем в общий сигнал выход алгоритма калибровки
    start_index_s = 1;
    end_index_s = sim_options.M*N;

    sig_adc = zeros(sim_options.M*N*j,1);
    for i = 1:j
        sig_adc(start_index_s:end_index_s) = x_tilda(:,i);
        start_index_s = end_index_s + 1;
        end_index_s = start_index_s+(sim_options.M*N)-1;
    end

    %% Графики 

    % График RMS delta
    figure(7);
    for i = 1:sim_options.M
        subplot(sim_options.M,1,i)
        plot(rms_delta(i,:));
        title(['RMS delta' num2str(i)])
        xlabel('Номер значения') 
        ylabel('Значение') 
    end

    % График корреляций
    figure(8);
    for i = 1:sim_options.M
        subplot(sim_options.M,1,i)
        plot(cor(i,:));
        title(['Корреляция суб-АЦП' num2str(i)])
        xlabel('Номер итерации') 
        ylabel('Значение') 
    end

    % Графики выходных сигналов
    figure(9)
    subplot(sim_options.M+4,1,1)
    plot([s_to_subadc_int(end-250:end), s_after_subadc(end-250:end), sig_adc(end-250:end)]);
    title('Входной сигнал и выходной сигнал')
    xlabel('Номер отсчета') 
    ylabel('Значение отсчета') 
    legend({'Вход без ошибок','Вход с ошибками', 'Выход'},'Location','northeast');
    % дельты
    for i = 1:sim_options.M
        subplot(sim_options.M+4,1,i+1)
        plot(delta_tilda(i,:));
        title(['Дельта суб-АЦП' num2str(i)])
        xlabel('Номер итерации') 
        ylabel('Значение дельты') 
    end

    % ЧХ сигналов
    subplot(sim_options.M+4,1,sim_options.M+2)
    snr(s_to_subadc_int*2^-11, sim_options.Fs_sub_adc*sim_options.M);
    subplot(sim_options.M+4,1,sim_options.M+3)
    snr(s_after_subadc*2^-11, sim_options.Fs_sub_adc*sim_options.M);
    subplot(sim_options.M+4,1,sim_options.M+4)
	snr(sig_adc*2^-11, sim_options.Fs_sub_adc*sim_options.M);

end