
% 1) Hu.M, Yi.P, (2022), Digital Calibration for Gain, Time Skew, and Bandwidth Mismatch 
%    in Under-Sampling Time-Interleaved System
% 2) Джиган В.И, Адаптивные фильтры
% 3) Айфичер Э, Джервис Б, Цифровая обработка сигналов. Практический подход
% 4) Behrouz Farhang-Boroujeny, Adaptive Filters Theory and Applications 

function [x_after_adc, x_after_adc_double, x_after_adc_int, fractional_mult, fractional_sum, fractional_width_total_mult, fractional_width_total_sum, y_fractional_outInt_abs_max, ...  
    hilbert_mult, hilbert_sum, hilbert_width_total_mult, hilbert_width_total_sum, ymi_HilbertInt_abs_max ...
    ... % Determinant
	... % умножители определителя 2x2 
	DetM_2x2_multiplier_total_abs_max, ...
	... % сумматоры определителя 2x2 
	Det2x2_sum_abs_max, ...
	... % умножители определителя 3х3
	Mult_DetM_3x3_array_max, ...
	... % разрядность умножителей определителя 3x3
	Mult_DetM_3x3_array_mult_total_width_max, ...
	... % пресумматоры определителя 3х3
	DetM_3x3_int_pre_sum_array_max, ...
	... % разрядность пресумматоров определителя 3x3
	DetM_3x3_int_pre_sum_width_total_max, ...
	... % сумматоры определителя 3х3
	DetM_3x3_int_sum_array_max, ...
	... % разрядность сумматоров определителя 3x3
	DetM_3x3_int_sum_width_total_max, ...
	... % умножители определителя 4х4
	DetM_4x4_int_mult_array_max, ...
	... % разрядность умножителей определителя 4х4
	DetM_4x4_int_mult_width_total_max, ...
	... % пресумматоры определителя 4х4
	DetM_4x4_int_pre_sum_array_max, ...
	... % разрядность пресумматоров определителя 4х4
	DetM_4x4_int_pre_sum_width_total_max, ...
	... % сумматоры определителя 4х4
	DetM_4x4_int_sum_array_max, ...
	... % разрядность сумматоров определителя 4х4
	DetM_4x4_int_sum_array_width_total_max, ...
	... % умножители определителя 5х5
	DetM_5x5_int_mult_array_max, ...
	... % разрядность умножителей 5x5
	DetM_5x5_int_mult_array_width_total_max, ...
	... % пресумматоры1 определителя 5х5
	DetM_5x5_int_pre_sum1_array_max, ...
	... % разрядность пресумматоров1 определителя 5х5
	DetM_5x5_int_pre_sum1_array_width_total_max, ...
	... % пресумматор2 определителя 5х5
	DetM_5x5_int_sum3_abs_max, ...
	... % разрядность пресумматора2 определителя 5х5
	DetM_5x5_int_sum3_width_total_max, ...
	... % сумматор определителя 5х5
	DetM_5x5_int_abs_max, ...
	... % разрядность сумматора определителя 5х5
	DetM_5x5_int_width_total_max, ...
    ... % начальный определитель
    Det_x3_int_max, ...
    ... % выход делителя
    Divide_max, ...
    ... % умножители адаптивного фильтра
    Adaptive_filter_mult_array_max, ...
    ... % разрядность умножителей адаптивного фильтра
    Adaptive_filter_mult_total_width, ...
    ... % сумматоры адаптивного фильтра
    Adaptive_filter_sum_array_max, ....
    ... % разрядность сумматоров адаптивного фильтра
    Adaptive_filter_sum_total_width ...
] = adc_calibration(sim_options, adc_input, s_to_subadc_int, s_after_subadc)

    %% Calibration algorithm 1 (Fractional delays)
    n = (0:1:sim_options.N-1);
    Nbp = floor(sim_options.Z/2); % стр 7. (24)
    del_proc = ((sim_options.N-1)/2);
    nn = 1:length(adc_input(:,1));
    w_blackman = 0.42 - 0.5 * cos(2*pi*n/(sim_options.N-1)) + 0.08 * cos(4*pi*n/(sim_options.N-1)); % Blackman window

    %% Fractional filter coeff
    %%
    delay_adc = (1/sim_options.M:1/sim_options.M:1); % (стр.6,(16)), создаем массив на различные значения задержек
    D = del_proc - delay_adc; % delay (N-1)/2 - d = causal filter
    hri_m = sinc(n'- D); % shift impulse response on D = Dint - d for fractional delay filter
    w_blackman_fractional = 0.42 - 0.5 * cos(2*pi*(n'+ delay_adc)/(sim_options.N-1)) + 0.08 * cos(4*pi*(n'+ delay_adc)/(sim_options.N-1)); % shift Blackman window
    hri_w = hri_m .* w_blackman_fractional; 
  
    coeff_frac_int = cast((hri_w*2^(sim_options.fractional_coeff_width-1)), sim_options.int_size);
    % запись коэффициентов фильтра дробной задержки в файл
    % for i = 1:sim_options.M-1
    %     writematrix(coeff_frac_int(:,i), ['src/width_txt/Коэффициенты_фильтров_дробной_задержки_АЦП' num2str(i), '.txt']);
    % end

    % figure(3)
    % subplot(2,1,1)
    % plot(double(coeff_frac_int))
    % title('Импульсная характеристика фильтра Гилберта')
    % xlabel('Номер отсчета') 
    % ylabel('Амплитуда') 
    % subplot(2,1,2)
    % plot(width)
    % title('Разрядность коэффициентов')
    % xlabel('Номер коэффициента') 
    % ylabel('Необходимое количество бит') 

    for i = 1:sim_options.M
        [y3(:,i), f3(:,i)] = freqz(hri_w(:,i), 1,1024, 'whole', 1000000000);
        [y4(:,i), f4(:,i)] = freqz(double(coeff_frac_int(:,i))*2^-(sim_options.fractional_coeff_width-1),1,1024, 'whole', 1000000000);
    end
    
    figure(2);
    subplot(4,1,1)
    plot(f3(:,1), abs(y3(:,1)), f3(:,1), abs(y4(:,1)));
    subplot(4,1,2)
    plot(f3(:,2), abs(y3(:,2)), f3(:,2), abs(y4(:,2)));
    subplot(4,1,3)
    plot(f3(:,3), abs(y3(:,3)), f3(:,3), abs(y4(:,3)));
    subplot(4,1,4)
    plot(f3(:,4), abs(y3(:,4)), f3(:,4), abs(y4(:,4)));

    % title('Влияние разрядностей коэффициентов на АЧХ фильтра дробной задержки')
    % xlabel('Частота') 
    % ylabel('Коэффициент передачи') 
    % legend({'double','16 бит', '19 бит', '21 бит'},'Location','northeast')

    %% Hilbert filter coeff
    %%
    hh = (2./((n-del_proc)*pi)).*(sin(((n-del_proc)*pi)./2)).^2;
    hh(1) = 0;
    hh(37) = 0;
    hh_m = (hh .* w_blackman).';

    % Negative Symmetric coefficients
    hilbert_coeff_int = cast(hh_m*2^(sim_options.hilbert_coeff_width-1), sim_options.int_size);
    % запись коэффициентов фильтра Гилберта в файл
    % writematrix(hilbert_coeff_int, ['src/width_txt/Коэффициенты_фильтра_Гилберта.txt']);

    % [y, f] = freqz(double(hilbert_coeff_int)*2^-(sim_options.hilbert_coeff_width-1), 1,1024, 'whole', 1000000000);
    % [y1, f1] = freqz(hh_m, 1,1024, 'whole', 1000000000);
    % figure(3);
    % plot(f, abs(y), f, abs(y1));
    
    %% zones Nyquist
    nn1 = nn' + delay_adc;
    a1 = 2*pi*nn1*Nbp;
    cosi = cos(a1);
    sini = sin(a1);

    %%
    yri_cut = zeros(length(adc_input(1:end-del_proc,1)),sim_options.M);
    yri_cut_int = cast(zeros(length(adc_input(1:end-del_proc,1)),sim_options.M), sim_options.int_size);

    yric = cast(zeros(length(adc_input(:,1)),sim_options.M-1), sim_options.int_size);
    yric_int = cast(zeros(length(adc_input(:,1)),sim_options.M-1), sim_options.int_size);

    fractional_mult = cast(zeros(sim_options.N, sim_options.M-1), sim_options.int_size); 
    fractional_sum = cast(zeros(sim_options.N-1, sim_options.M-1), sim_options.int_size); 
    fractional_width_total_mult = cast(zeros(sim_options.N, sim_options.M-1), sim_options.int_size); 
    fractional_width_total_sum = cast(zeros(sim_options.N-1, sim_options.M-1), sim_options.int_size); 

    hilbert_mult = cast(zeros(sim_options.N, sim_options.M-1), sim_options.int_size); 
    hilbert_sum = cast(zeros(sim_options.N-1, sim_options.M-1), sim_options.int_size); 
    hilbert_width_total_mult = cast(zeros(sim_options.N, sim_options.M-1), sim_options.int_size); 
    hilbert_width_total_sum = cast(zeros(sim_options.N-1, sim_options.M-1), sim_options.int_size); 

    y_fractional_outInt_div = cast(zeros(length(adc_input(:,1)),1), sim_options.type_fir_out);
    y_fractional_outInt_abs = cast(zeros(length(adc_input(:,1)),1), sim_options.type_fir_out);
    y_fractional_outInt_abs_max = cast(zeros(sim_options.M-1,1), sim_options.type_fir_out);

    ymi_HilbertInt_div = cast(zeros(length(adc_input(:,1)),1), sim_options.type_fir_out);
    ymi_HilbertInt_abs = cast(zeros(length(adc_input(:,1)),1), sim_options.type_fir_out);
    ymi_HilbertInt_abs_max = cast(zeros(sim_options.M-1,1), sim_options.type_fir_out);
    
	%% Задержка отсчетов сигнала АЦП0
    for i = 1:sim_options.M-1

        fractional_mult_file = ['src/width_txt/Width_multiplier_Fractional_filter_' num2str(i) '.txt'];
        fractional_sum_file = ['src/width_txt/Width_adder_Fractional_filter_' num2str(i) '.txt'];

        hilbert_mult_file = ['src/width_txt/Width_multiplier_Hilbert_filter_' num2str(i) '.txt'];
        hilbert_sum_file = ['src/width_txt/Width_adder_Hilbert_filter_' num2str(i) '.txt'];

        yri = filter(hri_w(:,i), 1, adc_input(:,1)); % filter (стр.6 (15))
        yri = [yri(del_proc+1:end); zeros(del_proc,1)]; % убираем переходной процесс

        yri_round = round(yri);
        % yri_floor = floor(yri);
        % yri_ceil = ceil(yri);
        % yri_fix = fix(yri);

        % figure(2);
        % subplot(5,1,1)
        % snr(yri, 1000000000);
        % subplot(5,1,2);
        % snr(yri_round, 1000000000);
        % subplot(5,1,3);
        % snr(yri_floor, 1000000000);
        % subplot(5,1,4);
        % snr(yri_ceil, 1000000000);
        % subplot(5,1,5);
        % snr(yri_fix, 1000000000);

        ymi = filter(hh_m.', 1, yri);
        ymi = [ymi(del_proc+1:end); zeros(del_proc,1)]; % убираем переходной процесс

        ymi_round = round(ymi);

        %% Фильтр дробной задержки
        [y_fractional_outInt, fractional_mult(:,i), fractional_sum(:,i), fractional_width_total_mult(:,i), fractional_width_total_sum(:,i) ...
            ] = fir_filter(coeff_frac_int(:,i), adc_input(:,1), sim_options.N, ...
        fractional_mult_file, fractional_sum_file, sim_options.width_fractional, sim_options); % (стр.6 (15)) 

        y_fractional_outInt = [y_fractional_outInt(del_proc+1:end); zeros(del_proc,1)]; % убираем переходной процесс

        for j = 1:length(y_fractional_outInt)
            % округляем инты
            if (bitget(y_fractional_outInt(j), (sim_options.fractional_coeff_width-1)) == 1)
                y_fractional_outInt_div(j) = bitshift(y_fractional_outInt(j), -(sim_options.fractional_coeff_width-1)); % сдвигаем данные
                y_fractional_outInt_div(j) = y_fractional_outInt_div(j) + cast(1,sim_options.type_fir_out);
            else
                y_fractional_outInt_div(j) = bitshift(y_fractional_outInt(j), -(sim_options.fractional_coeff_width-1)); % сдвигаем данные
            end

            % возвращаем поделенные положительные значения 
            if y_fractional_outInt_div(j) < 0
                y_fractional_outInt_abs(j) = y_fractional_outInt_div(j) * cast(-1, sim_options.type_fir_out); % находим число по модулю
            else
	            y_fractional_outInt_abs(j) = y_fractional_outInt_div(j);
            end

            % ищем максимум
            if (y_fractional_outInt_abs_max(i,1) < y_fractional_outInt_abs(j))
                y_fractional_outInt_abs_max(i,1) = y_fractional_outInt_abs(j);
            end
        end

        y_fractional_outInt_double = double(y_fractional_outInt_div);
        relative_error_fractional = yri_round./y_fractional_outInt_double;

        figure(4);
        subplot(5,1,1)
        plot([yri_round(1:900), y_fractional_outInt_double(1:900)]);
        subplot(5,1,2);
        snr(yri, 1000000000);
        subplot(5,1,3);
        snr(yri_round, 1000000000);
        subplot(5,1,4);
        snr(y_fractional_outInt_double, 1000000000);
        subplot(5,1,5);
        plot(relative_error_fractional);
        title('Относительная ошибка выходного сигнала фильтра дробной задержки между double и integer')
        xlabel('Номер отсчета') 
        ylabel('Значение ошибки') 
        x4 = xline(37, '--', 'Переходной процесс фильтра')
        x4.LabelHorizontalAlignment = 'center'
        x4.LabelVerticalAlignment = 'middle';

        %% Фильтр Гилберта
        
        [ymi_HilbertInt, hilbert_mult(:,i), hilbert_sum(:,i), hilbert_width_total_mult(:,i), hilbert_width_total_sum(:,i)] = ...
            fir_filter(hilbert_coeff_int, y_fractional_outInt_div, sim_options.N, hilbert_mult_file, hilbert_sum_file, sim_options.width_hilbert, sim_options); % (стр.6 (15)) );
        ymi_HilbertInt = [ymi_HilbertInt(del_proc+1:end); zeros(del_proc,1)]; % убираем переходной процесс
 
        % округляем инты
        for j = 1:length(ymi_HilbertInt)
            if (bitget(ymi_HilbertInt(j), (sim_options.hilbert_coeff_width-1)) == 1)
                ymi_HilbertInt_div(j) = bitshift(ymi_HilbertInt(j), -(sim_options.hilbert_coeff_width-1)); % сдвигаем данные
                ymi_HilbertInt_div(j) = ymi_HilbertInt_div(j) + cast(1,sim_options.type_fir_out);
            else
                ymi_HilbertInt_div(j) = bitshift(ymi_HilbertInt(j), -(sim_options.hilbert_coeff_width-1)); % сдвигаем данные
            end

            % возвращаем поделенные положительные значения 
            if ymi_HilbertInt_div(j) < 0
                ymi_HilbertInt_abs(j) = ymi_HilbertInt_div(j) * cast(-1, sim_options.type_fir_out); % находим число по модулю
            else
	            ymi_HilbertInt_abs(j) = ymi_HilbertInt_div(j);
            end

            % ищем максимум
            if (ymi_HilbertInt_abs_max(i,1) < ymi_HilbertInt_abs(j))
                ymi_HilbertInt_abs_max(i,1) = ymi_HilbertInt_abs(j);
            end
        end

        ymi_HilbertInt_double = double(ymi_HilbertInt_div);
        relative_error_hilbert = ymi_round./ymi_HilbertInt_double;

        figure(5);
        subplot(5,1,1)
        plot([ymi_round(1:500), double(ymi_HilbertInt_div(1:500))]);
        subplot(5,1,2);
        snr(ymi, 1000000000);
        subplot(5,1,3);
        snr(ymi_round, 1000000000);
        subplot(5,1,4);
        snr(double(ymi_HilbertInt_div), 1000000000);
        subplot(5,1,5);
        plot(relative_error_hilbert);
        title('Относительная ошибка выходного сигнала фильтра Гилберта между double и integer')
        xlabel('Номер отсчета') 
        ylabel('Значение ошибки') 
        x4 = xline(73, '--', 'Переходной процесс фильтра')
        x4.LabelHorizontalAlignment = 'center'
        x4.LabelVerticalAlignment = 'middle';

        %%
        [yric(:,i), yric_int(:,i)] = single_sideband(yri, y_fractional_outInt_div, ymi, ymi_HilbertInt_div, cosi(:,i), sini(:,i), i, sim_options);

    end

    for i = 1:sim_options.M
        if i == 1
             yri_cut(:,1) = adc_input(1:end-del_proc,1);
             yri_cut_int(:,1) = adc_input(1:end-del_proc,1);
        else
            yri_cut(:,i) = yric(1:end-del_proc,i-1);
            yri_cut_int(:,i) = yric_int(1:end-del_proc,i-1);
        end
    end

    % test signal
    sig_adc = zeros(sim_options.M*length(yri_cut(:,1)),1);
    sig_adc_int = zeros(sim_options.M*length(yri_cut_int(:,1)),1);
	for i = 1:sim_options.M
        sig_adc(i:sim_options.M:end) = yri_cut(:,i);
        sig_adc_int(i:sim_options.M:end) = double(yri_cut_int(:,i));
    end

    %%
    % save (sprintf(num2str(clock) + ".mat"));
    % load ('2025             11             20             12             29         36.135.mat'); 
    % unique
    % load ('2025             11             13             17              3           49.5.mat'); % 50 MHz 70 SNR

    figure(15);
    subplot(2,1,1)
    plot([sig_adc(1:750), sig_adc_int(1:750)]);
    subplot(2,1,2)
    plot(sig_adc ./ sig_adc_int);

   figure(16);
   subplot(3,1,1)
   snr(s_to_subadc_int, sim_options.Fs/sim_options.Inter);
   subplot(3,1,2);
   snr(sig_adc, sim_options.Fs/sim_options.Inter);
   subplot(3,1,3);
   snr(sig_adc_int, sim_options.Fs/sim_options.Inter);

   % sim_options.type_2x2_det =      "int64";                ... % тип данных для матрицы 2х2
   % sim_options.type_3x3_det =      "int64";                ... % тип данных для матрицы 3х3
   % sim_options.type_4x4_det =      "int64";                ...    
   % sim_options.type_5x5_det =      "int64";                ...
   % sim_options.type_divide_out =  "int64";
   % sim_options.divide_factor = 14;
   % sim_options.type_mult_in_adaptive_filter    = "int64";
   % sim_options.type_add_in_adaptive_filter     = "int64";
   % 
   % sim_options.type_2x2_det =      "int64";                ... % тип данных для матрицы 2х2
   %% Calibration algorithm 2 (Least Mean Squares)
   for i = 1:sim_options.M-1
        adaptive_mult_file = ['src/width_txt/Width_mult_Adaptive_filter_' num2str(i) '.txt'];
        adaptive_sum_file = ['src/width_txt/Width_adder_Adaptive_filter_' num2str(i) '.txt'];

        [y_array(:,i), y_array_double(:,i), y_array_int(:,i), DetM_2x2_array(:,i), DetM_2x2_array_int(:,i), ...
			... % умножители определителя 2x2 
			DetM_2x2_multiplier_total_abs_max(:,i), ...
			... % сумматоры определителя 2x2 
			Det2x2_sum_abs_max(:,i), ...
			... % умножители определителя 3х3
			Mult_DetM_3x3_array_max(:,i), ...
			... % разрядность умножителей определителя 3x3
			Mult_DetM_3x3_array_mult_total_width_max(:,i), ...
			... % пресумматоры определителя 3х3
			DetM_3x3_int_pre_sum_array_max(:,i), ...
			... % разрядность пресумматоров определителя 3x3
			DetM_3x3_int_pre_sum_width_total_max(:,i), ...
			... % сумматоры определителя 3х3
			DetM_3x3_int_sum_array_max(:,i), ...
			... % разрядность сумматоров определителя 3x3
			DetM_3x3_int_sum_width_total_max(:,i), ...
			... % умножители определителя 4х4
			DetM_4x4_int_mult_array_max(:,i), ...
			... % разрядность умножителей определителя 4х4
			DetM_4x4_int_mult_width_total_max(:,i), ...
			... % пресумматоры определителя 4х4
			DetM_4x4_int_pre_sum_array_max(:,i), ...
			... % разрядность пресумматоров определителя 4х4
			DetM_4x4_int_pre_sum_width_total_max(:,i), ...
			... % сумматоры определителя 4х4
			DetM_4x4_int_sum_array_max(:,i), ...
			... % разрядность сумматоров определителя 4х4
			DetM_4x4_int_sum_array_width_total_max(:,i), ...
			... % умножители определителя 5х5
			DetM_5x5_int_mult_array_max(:,i), ...
			... % разрядность умножителей 5x5
			DetM_5x5_int_mult_array_width_total_max(:,i), ...
			... % пресумматоры1 определителя 5х5
			DetM_5x5_int_pre_sum1_array_max(:,i), ...
			... % разрядность пресумматоров1 определителя 5х5
			DetM_5x5_int_pre_sum1_array_width_total_max(:,i), ...
			... % пресумматор2 определителя 5х5
			DetM_5x5_int_sum3_abs_max(:,i), ...
			... % разрядность пресумматора2 определителя 5х5
			DetM_5x5_int_sum3_width_total_max(:,i), ...
			... % сумматор определителя 5х5
			DetM_5x5_int_abs_max(:,i), ...
			... % разрядность сумматора определителя 5х5
			DetM_5x5_int_width_total_max(:,i), ...
            ... % начальный определитель
            Det_x3_int_max(:,i), ...
            ... % выход делителя
            Divide_max(:,i), ...
            ... % умножители адаптивного фильтра
            Adaptive_filter_mult_array_max(:,i), ...
            ... % разрядность умножителей адаптивного фильтра
            Adaptive_filter_mult_total_width(:,i), ...
            ... % сумматоры адаптивного фильтра
            Adaptive_filter_sum_array_max(:,i), ....
            ... % разрядность сумматоров адаптивного фильтра
            Adaptive_filter_sum_total_width(:,i) ...
        ] = least_mean_squares(adc_input(:,i+1), yri_cut(:,i+1), yri_cut_int(:,i+1), adaptive_mult_file, adaptive_sum_file, sim_options);
    end

    % save (sprintf(num2str(clock) + ".mat"));
    % load ('2025             11             12             12             40         40.919.mat'); % 888 MHz 70 SNR

    %% create main signal after LS algorithm (switch after sub-adc)
    x_after_adc = zeros(length(y_array)*sim_options.M,1);
    x_after_adc_double = zeros(length(y_array_double)*sim_options.M,1);
    x_after_adc_int = zeros(length(y_array_int)*sim_options.M,1);
    
    for i = 1:sim_options.M
        if i == 1
            x_after_adc(i:sim_options.M:end) = yri_cut(1:length(y_array),1);
            x_after_adc_double(i:sim_options.M:end) = yri_cut(1:length(y_array_double),1);
            x_after_adc_int(i:sim_options.M:end) = double(yri_cut_int(1:length(y_array_int),1));
        else
            x_after_adc(i:sim_options.M:end) = y_array(:,i-1);
            x_after_adc_double(i:sim_options.M:end) = y_array_double(:,i-1); % * 2^-(sim_options.divide_factor);
            x_after_adc_int(i:sim_options.M:end) = double(y_array_int(:,i-1)); % * 2^-(sim_options.divide_factor);
        end
    end

    figure(7);
    subplot(3,1,1)
    plot([s_to_subadc_int(1:500), x_after_adc(1:500)]);
    title('Исходный сигнал и выход алгоритма LU')
    xlabel('Номер отсчета') 
    ylabel('Амплитуда') 
    legend({'Исходный сигнал','Сигнал с выхода адаптивного фильтра double'},'Location','northeast')
    subplot(3,1,2)
    plot([s_to_subadc_int(1:600), x_after_adc_double(1:600)]);
    title('Исходный сигнал и выход адаптивного фильтра double')
    xlabel('Номер отсчета') 
    ylabel('Амплитуда') 
    legend({'Исходный сигнал','Сигнал с выхода адаптивного фильтра '},'Location','northeast')
    subplot(3,1,3)
    plot([s_to_subadc_int(1:500), x_after_adc_int(1:500)]);
    title('Исходный сигнал и выход адаптивного фильтра int')
    xlabel('Номер отсчета') 
    ylabel('Амплитуда');
    legend({'Исходный сигнал','Сигнал с выхода адаптивного фильтра int'},'Location','northeast')

end
