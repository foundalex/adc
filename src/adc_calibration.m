
% 1) Hu.M, Yi.P, (2022), Digital Calibration for Gain, Time Skew, and Bandwidth Mismatch 
%    in Under-Sampling Time-Interleaved System
% 2) Джиган В.И, Адаптивные фильтры
% 3) Айфичер Э, Джервис Б, Цифровая обработка сигналов. Практический подход
% 4) Behrouz Farhang-Boroujeny, Adaptive Filters Theory and Applications 

function [x_after_adc, x_after_adc_double, x_after_adc_int, ...
    ... % Значения полосового фильтра ADC0
    filter_max_width_out_golden, ...
    ... % Значения фильтров дробной задержки
    filter_max_width_out_fractional, ...
    y_fractional_outInt_abs_max, ...
    ... % Значения фильтров Гилберта
    filter_max_width_out_hilbert, ... % hilbert_mult, hilbert_sum, hilbert_width_total_mult, hilbert_width_total_sum, ymi_HilbertInt_abs_max, ...
    ymi_HilbertInt_abs_max, ...
    ... % Determinant
    determinate_struct, ...
    ... % выход делителя
    Divide_max, ...
    adaptive_filter_struct_max ...
] = adc_calibration(sim_options, adc_input, s_to_subadc_int, s_after_subadc)

    %% Различные переменные
    % Переменные полосового фильтра АЦП0
    y_golden_outInt_abs = cast(zeros(length(adc_input(:,1)),1), sim_options.type_fir_out);
    
    % Переменные фильтров дробной задержки
    y_fractional_outInt_abs = cast(zeros(length(adc_input(:,1)),1), sim_options.type_fir_out);
    y_fractional_outInt_abs_max = cast(zeros(sim_options.M-1,1), sim_options.type_fir_out);

    % Переменные фильтра Гилберта
    ymi_HilbertInt_abs = cast(zeros(length(adc_input(:,1)),1), sim_options.type_fir_out);
    ymi_HilbertInt_abs_max = cast(zeros(sim_options.M-1,1), sim_options.type_fir_out);

    % Переменные модуля переноса частоты
    yric = cast(zeros(length(adc_input(:,1)),sim_options.M-1), sim_options.int_size);
    yric_int = cast(zeros(length(adc_input(:,1)),sim_options.M-1), sim_options.int_size);
  
    %%                                                       
    n = (0:1:sim_options.N-1);
    Nbp = floor(sim_options.Z/2); % стр 7. (24)
    del_proc = ((sim_options.N-1)/2);
    nn = 1:length(adc_input(:,1));
    % Окно Блэкмена
    w_blackman = 0.42 - 0.5 * cos(2*pi*n/(sim_options.N-1)) + 0.08 * cos(4*pi*n/(sim_options.N-1));

    %% Расчет коэффициентов полосовых фильтров дробной задержки
    %%
    % Частота среза 1-го ФНЧ
    ws1 = 0.91;
    % Частота среза 2-го ФНЧ
    ws2 = 0.04;
    % Массив различных значений задержек для фильтров дробной задержки
    delay_adc = (1/sim_options.M:1/sim_options.M:1); % (стр.6,(16))

    D = del_proc - delay_adc; % delay (N-1)/2 - d = causal filter
    % Синтезируем ФНЧ фильтры дробной задержки, из которых в дальнейшем
    % получим полосовые
    lowpass1_fractional = ws1*sinc(ws1*(n'- D)); % shift impulse response on D = Dint - d for fractional delay filter
    lowpass2_fractional = ws2*sinc(ws2*(n'- D)); % shift impulse response on D = Dint - d for fractional delay filter
    % Умножаем на сдвинутое окно Блэкмена
    w_blackman_fractional = 0.42 - 0.5 * cos(2*pi*(n'+ delay_adc)/(sim_options.N-1)) + 0.08 * cos(4*pi*(n'+ delay_adc)/(sim_options.N-1)); 
    % Из двух фильтров ФНЧ с разными частотами среза получаем полосовой,
    % путем вычитания одного из другого
    weight_lowpass1_fractional = lowpass1_fractional .* w_blackman_fractional; 
    weight_lowpass2_fractional = lowpass2_fractional .* w_blackman_fractional; 
    bandpass_fractional = weight_lowpass1_fractional - weight_lowpass2_fractional;

    % перевод коэффициентов фильтров дробной задержки в инты
    coeff_frac_int = cast((bandpass_fractional*2^(sim_options.fractional_coeff_width-1)), sim_options.int_size);

    % запись коэффициентов фильтра дробной задержки в файл 
    % for i = 1:sim_options.M-1
    %     writematrix(coeff_frac_int(:,i), ['src/width_txt/Коэффициенты_фильтров_дробной_задержки_АЦП' num2str(i), '.txt']);
    % end
    %% Синтезируем полосовой фильтр для АЦП0
    lowpass1_adc0 = ws1*sinc(ws1*(n'-del_proc));
    lowpass2_adc0 = ws2*sinc(ws2*(n'-del_proc));
   
    weight_lowpass1_adc0 = lowpass1_adc0 .* w_blackman'; 
    weight_lowpass2_adc0 = lowpass2_adc0 .* w_blackman'; 
    bandpass_adc0 = weight_lowpass1_adc0  - weight_lowpass2_adc0;

    % перевод коэффициентов эталонного фильтра в инты
    coeff_gold_adc0_int = cast((bandpass_adc0*2^(sim_options.hilbert_coeff_width-1)), sim_options.int_size);

    % [badc0_double_value, badc0_double_frequency] = freqz(bandpass_adc0, 1,1024, 'whole', sim_options.Fs_sub_adc);
    % [badc0_int_value, badc0_int_frequency] = freqz(double(coeff_gold_adc0_int), 1,1024, 'whole', sim_options.Fs_sub_adc);

    % Сравнение коэффициентов double и int
    % figure(2);
    % plot(badc0_double_frequency, abs(badc0_double_value), badc0_double_frequency, abs(badc0_int_value*(2^-(sim_options.hilbert_coeff_width-1))));
    % title('АЧХ полосового эталонного фильтра')
    % xlabel('Частота') 
    % ylabel('Коэффициент передачи') 
    % legend({'double', 'int'},'Location','northeast');
    % x1 = xline(500000000, '--', 'Fs/2')
    % x1.LabelHorizontalAlignment = 'center'
    % x1.LabelVerticalAlignment = 'middle';

    %% Частотная характеристика полосовых фильтров 
    % for i = 1:sim_options.M-1
    %     [y3(:,i), f3(:,i)] = freqz(bandpass_fractional(:,i), 1,1024, 'whole', sim_options.Fs_sub_adc);
    %     % [y4(:,i), f4(:,i)] = freqz(double(coeff_frac_int(:,i))*2^-(sim_options.fractional_coeff_width-1),1,1024, 'whole', sim_options.Fs_sub_adc);
    % end

    % figure(3);
    % plot(f3(:,1), abs(y3(:,1)), f3(:,1), abs(y3(:,2)), f3(:,1), abs(y3(:,3)), f3(:,1), abs(badc0_double_value));
    % title('АЧХ полосовых фильтров')
    % xlabel('Частота') 
    % ylabel('Коэффициент передачи') 
    % legend({'0.25*Fs','0.5*Fs', '0.75*Fs', 'Эталонный фильтр'},'Location','northeast');

    % x2 = xline(500000000, '--', 'Fs/2')
    % x2.LabelHorizontalAlignment = 'center'
    % x2.LabelVerticalAlignment = 'middle';

    %% Расчет коэффициентов фильтра Гилберта

    lowpass_hilbert1 = (2./((n-del_proc)*pi)).*(sin(ws1*((n-del_proc)*pi)./2)).^2;
    lowpass_hilbert1(1) = 0;
    lowpass_hilbert1(37) = 0;
    lowpass_hilbert2 = (2./((n-del_proc)*pi)).*(sin(ws2*((n-del_proc)*pi)./2)).^2;
    lowpass_hilbert2(1) = 0;
    lowpass_hilbert2(37) = 0;

    weight_lowpass_hilbert1 = lowpass_hilbert1 .* w_blackman; 
    weight_lowpass_hilbert2 = lowpass_hilbert2 .* w_blackman; 
    bandpass_hilbert = weight_lowpass_hilbert1'  - weight_lowpass_hilbert2';

    % Перевод коэффициентов фильтра Гилберта в инты
    hilbert_coeff_int = cast(bandpass_hilbert*2^(sim_options.hilbert_coeff_width-1), sim_options.int_size);
    % запись коэффициентов фильтра Гилберта в файл
    % writematrix(hilbert_coeff_int, ['src/width_txt/Коэффициенты_фильтра_Гилберта.txt']);

    % [y, f] = freqz(double(hilbert_coeff_int)*2^-(sim_options.hilbert_coeff_width-1), 1,1024, 'whole', 1000000000);
    % [y1, f1] = freqz(bandpass_hilbert, 1,1024, 'whole', 1000000000);

    % figure(4);
    % plot(f3(:,1), abs(y3(:,1)), f3(:,1), abs(y3(:,2)), f3(:,1), abs(y3(:,3)), f3(:,1), abs(badc0_double_value), f1, abs(y1));
    % title('АЧХ полосовых фильтров дробной задержки и фильтра Гилберта')
    % xlabel('Частота') 
    % ylabel('Коэффициент передачи') 
   
    %% Коэффициенты эталонного фильтра
    % golden_h = (fir1(128,[0.030 0.440],"bandpass"))';

    %% Данные для переноса сигнала в различные зоны Найквиста
    nn1 = nn' + delay_adc;
    a1 = 2*pi*nn1*Nbp;
    cosi = cos(a1);
    sini = sin(a1);

                                                                %% Первая часть алгоритма калибровки
    %% Эталонный сигнал
    % filter_max_width_out_golden = struct;

    % double
    yr_double = filter(bandpass_adc0, 1, adc_input(:,1));
    yr_double = [yr_double(del_proc+1:end); zeros(del_proc,1)]; 

    yr_double_round = round(yr_double);

    %% Int
    golden_max_width = readtable(['src/width_txt/Разрядность_максимальных_значений_эталонного_фильтра.xlsx']); 
    golden_max_width = table2array(golden_max_width(:,2:end));

    % Фильтруем эталонный сигнал
    [y_golden_outInt, filter_max_width_out_golden ...
        ] = fir_filter(coeff_gold_adc0_int, adc_input(:,1), sim_options.N, golden_max_width, sim_options.width_hilbert, sim_options);

    filter_max_width_out_golden.y_golden_outInt_abs_max = cast(zeros(1,1), sim_options.type_fir_out);

    % убираем переходной процесс
    y_golden_outInt = [y_golden_outInt(del_proc+1:end); zeros(del_proc,1)];

    % округляем значения
    y_golden_outInt_div = round_int(y_golden_outInt, sim_options.hilbert_coeff_width-1, sim_options.type_fir_out);

    for j = 1:length(y_golden_outInt_div)
        % возвращаем поделенные положительные значения 
        if y_golden_outInt_div(j) < 0
            y_golden_outInt_abs(j) = y_golden_outInt_div(j) * cast(-1, sim_options.type_fir_out); % находим число по модулю
        else
	        y_golden_outInt_abs(j) = y_golden_outInt_div(j);
        end

        % ищем максимум
        if (filter_max_width_out_golden.y_golden_outInt_abs_max < y_golden_outInt_abs(j))
            filter_max_width_out_golden.y_golden_outInt_abs_max = y_golden_outInt_abs(j);
        end
    end

    figure(5)
    subplot(5,1,1)
    plot([yr_double_round(1:900), y_golden_outInt_div(1:900)]);
    title('Выходные сигналы фильтра эталонного сигнала')
    xlabel('Номер отсчета') 
    ylabel('Значение отсчета') 
    subplot(5,1,2)
    snr((yr_double), sim_options.Fs_sub_adc);
    subplot(5,1,3)
    snr(yr_double_round, sim_options.Fs_sub_adc);
    subplot(5,1,4)
    snr(double(y_golden_outInt), sim_options.Fs_sub_adc);
    subplot(5,1,5)
    snr(double(y_golden_outInt_div), sim_options.Fs_sub_adc);

	%% Дробная задержка отсчетов сигнала с выхода АЦП0
    for i = 1:sim_options.M-1

        % fractional_filter_struct = struct;
        % filter_max_width_out_fractional = struct;

        fractional_max_width = readtable(['src/width_txt/Разрядность_максимальных_значений_фильтра_дробной_задержки_АЦП_№' num2str(i) '.xlsx']); 
        fractional_max_width = table2array(fractional_max_width(:,2:end));

        % fractional_mult_file = ['src/width_txt/Width_multiplier_Fractional_filter_' num2str(i) '.txt'];
        % fractional_sum_file = ['src/width_txt/Width_adder_Fractional_filter_' num2str(i) '.txt'];
        % 
        % hilbert_mult_file = ['src/width_txt/Width_multiplier_Hilbert_filter_' num2str(i) '.txt'];
        % hilbert_sum_file = ['src/width_txt/Width_adder_Hilbert_filter_' num2str(i) '.txt'];

        %% Фильтр дробной задержки (double)
        yri = filter(bandpass_fractional(:,i), 1, adc_input(:,1)); % filter (стр.6 (15))
        % убираем переходной процесс
        yri = [yri(del_proc+1:end); zeros(del_proc,1)]; % убираем переходной процесс
        % округляем значения
        yri_round = round(yri);

        %% Фильтр дробной задержки (Int)
       [y_fractional_outInt, filter_max_width_out_fractional(i)] = ...
            fir_filter(coeff_frac_int(:,i), adc_input(:,1), sim_options.N, fractional_max_width, sim_options.width_fractional, sim_options); % (стр.6 (15)) 

        % убираем переходной процесс
        y_fractional_outInt = [y_fractional_outInt(del_proc+1:end); zeros(del_proc,1)]; 
        % округляем значения
        y_fractional_outInt_div = round_int(y_fractional_outInt, sim_options.fractional_coeff_width-1, sim_options.type_fir_out);

        %% Ищем макс. значения для записи в файл
        for j = 1:length(y_fractional_outInt)
            % возвращаем поделенные положительные значения 
            if y_fractional_outInt_div(j) < 0
                y_fractional_outInt_abs(j) = y_fractional_outInt_div(j) * cast(-1, sim_options.type_fir_out); % находим число по модулю
            else
	            y_fractional_outInt_abs(j) = y_fractional_outInt_div(j);
            end

            % ищем максимум
            if (y_fractional_outInt_abs_max(i) < y_fractional_outInt_abs(j))
                y_fractional_outInt_abs_max(i) = y_fractional_outInt_abs(j);
            end
        end

        %% 
        y_fractional_outInt_double = double(y_fractional_outInt_div);
        relative_error_fractional = yri_round./y_fractional_outInt_double;

        figure(6);
        subplot(5,1,1)
        plot([yri_round(1:900), y_fractional_outInt_double(1:900)]);
        title('Выходные сигналы фильтров дробной задержки')
        xlabel('Номер отсчета') 
        ylabel('Значение отсчета') 
        legend({'double','int'},'Location','northeast')
        subplot(5,1,2);
        snr(yri, sim_options.Fs_sub_adc);
        subplot(5,1,3);
        snr(yri_round, sim_options.Fs_sub_adc);
        subplot(5,1,4);
        snr(y_fractional_outInt_double, sim_options.Fs_sub_adc);
        subplot(5,1,5);
        plot(relative_error_fractional);
        title('Относительная ошибка выходного сигнала фильтра дробной задержки между double и integer')
        xlabel('Номер отсчета') 
        ylabel('Значение ошибки') 
        x4 = xline(37, '--', 'Переходной процесс фильтра')
        x4.LabelHorizontalAlignment = 'center'
        x4.LabelVerticalAlignment = 'middle';

        %%
        % hilbert_filter_struct = struct;

        hilbert_max_width = readtable(['src/width_txt/Разрядность_максимальных_значений_фильтра_Гилберта_АЦП_№' num2str(i) '.xlsx']); 
        hilbert_max_width = table2array(hilbert_max_width(:,2:end));

        %% Фильтр Гилберта (Double)
        ymi = filter(bandpass_hilbert.', 1, yri);
        % убираем переходной процесс
        ymi = [ymi(del_proc+1:end); zeros(del_proc,1)]; 
        % округляем значения
        ymi_round = round(ymi);
        %% Фильтр Гилберта (Int)

        [ymi_HilbertInt, filter_max_width_out_hilbert(i)] = ...
            fir_filter(hilbert_coeff_int, y_fractional_outInt_div, sim_options.N, hilbert_max_width, sim_options.width_hilbert, sim_options); % (стр.6 (15)) );
        ymi_HilbertInt = [ymi_HilbertInt(del_proc+1:end); zeros(del_proc,1)]; % убираем переходной процесс

        % округляем значения
        ymi_HilbertInt_div = round_int(ymi_HilbertInt, sim_options.hilbert_coeff_width-1, sim_options.type_fir_out);

        %% Ищем макс. значения для записи в файл
        for j = 1:length(ymi_HilbertInt)
            % возвращаем поделенные положительные значения 
            if ymi_HilbertInt_div(j) < 0
                ymi_HilbertInt_abs(j) = ymi_HilbertInt_div(j) * cast(-1, sim_options.type_fir_out); % находим число по модулю
            else
	            ymi_HilbertInt_abs(j) = ymi_HilbertInt_div(j);
            end

            % ищем максимум
            if (ymi_HilbertInt_abs_max(i) < ymi_HilbertInt_abs(j))
                ymi_HilbertInt_abs_max(i) = ymi_HilbertInt_abs(j);
            end
        end

        %%
        ymi_HilbertInt_double = double(ymi_HilbertInt_div);
        relative_error_hilbert = ymi_round./ymi_HilbertInt_double;

        figure(7);
        subplot(5,1,1)
        plot([ymi_round(1:500), double(ymi_HilbertInt_div(1:500))]);
        title('Выходные сигналы фильтра Гилберта')
        xlabel('Номер отсчета') 
        ylabel('Значение отсчета') 
        legend({'double','int'},'Location','northeast')
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

        %% Перенос сигнала в заданные зоны Найквиста
        [yric(:,i), yric_int(:,i)] = single_sideband(yri, y_fractional_outInt_div, ymi, ymi_HilbertInt_div, cosi(:,i), sini(:,i), i, sim_options);

    end

    yri_cut = zeros(length(adc_input(1:end-del_proc,1)),sim_options.M);
    yri_cut_int = cast(zeros(length(adc_input(1:end-del_proc,1)),sim_options.M), sim_options.int_size);

    % Собираем полученные сигналы в массив для удобства
    for i = 1:sim_options.M
        if i == 1
             yri_cut(:,1) = yr_double(1:end-del_proc,1);
             yri_cut_int(:,1) = y_golden_outInt_div(1:end-del_proc,1);
        else
            yri_cut(:,i) = yric(1:end-del_proc,i-1);
            yri_cut_int(:,i) = yric_int(1:end-del_proc,i-1);
        end
    end

    %% Тестируем первую часть алгоритма калибровки
    sig_adc_gold = zeros(sim_options.M*length(adc_input(:,1)),1);
    sig_adc = zeros(sim_options.M*length(yri_cut(:,1)),1);
    sig_adc_int = zeros(sim_options.M*length(yri_cut_int(:,1)),1);
	for i = 1:sim_options.M
        sig_adc_gold(i:sim_options.M:end) = adc_input(:,i);
        sig_adc(i:sim_options.M:end) = yri_cut(:,i);
        sig_adc_int(i:sim_options.M:end) = double(yri_cut_int(:,i));
    end
    % 
    % save (sprintf(num2str(clock) + ".mat"));
    % load ('2025             11             25             17             39          22.88.mat'); 

    figure(8);
    subplot(3,1,1)
    plot([sig_adc_gold(1:750), sig_adc(1:750)]);
    subplot(3,1,2)
    plot([sig_adc(1:750), sig_adc_int(1:750)]);
    subplot(3,1,3)
    plot(sig_adc ./ sig_adc_int);

    figure(9);
    subplot(4,1,1)
    snr(s_to_subadc_int, sim_options.Fs/sim_options.Inter);
    subplot(4,1,2)
    snr(s_after_subadc, sim_options.Fs/sim_options.Inter);
    subplot(4,1,3);
    snr(sig_adc, sim_options.Fs/sim_options.Inter);
    subplot(4,1,4);
    snr(sig_adc_int, sim_options.Fs/sim_options.Inter);

    figure(10);
    subplot(4,1,1)
    sfdr(s_to_subadc_int, sim_options.Fs/sim_options.Inter);
    subplot(4,1,2)
    sfdr(s_after_subadc, sim_options.Fs/sim_options.Inter);
    subplot(4,1,3);
    sfdr(sig_adc, sim_options.Fs/sim_options.Inter);
    subplot(4,1,4);
    sfdr(sig_adc_int, sim_options.Fs/sim_options.Inter);


                                                   %% Вторая часть алгоритма калибровки (Метод наименьших квадратов)
   read_max_width = struct;

   for i = 1:sim_options.M-1

        det_max_width = readtable(['src/width_txt/Разрядность_максимальных_значений_определителя_АЦП_№' num2str(i) '.xlsx']); 
        read_max_width.det_max_width = table2array(det_max_width(:,2:end));

        adaptive_max_width = readtable(['src/width_txt/Разрядность_максимальных_значений_адаптивного_фильтра_АЦП_№' num2str(i) '.xlsx']); 
        read_max_width.adaptive_max_width = table2array(adaptive_max_width(:,2:end));

        [y_array(:,i), y_array_double(:,i), y_array_int(:,i), DetM_2x2_array(:,i), DetM_2x2_array_int(:,i), ...
			determinate_struct(i), ...
            Divide_max(i), ...
            adaptive_filter_struct_max(i) ...
        ] = least_mean_squares(adc_input(:,i+1), yri_cut(:,i+1), yri_cut_int(:,i+1), read_max_width, sim_options);
        
   end

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

    figure(11);
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
