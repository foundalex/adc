
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
] = adc_calibration(sim_options, adc_input_double, adc_input, s_to_subadc_int, s_after_subadc)

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
    yric = zeros(length(adc_input(:,1)),sim_options.M-1);
    yric_int = cast(zeros(length(adc_input(:,1)),sim_options.M-1), sim_options.type_fir_out);
  
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
    ws1 = 0.998;
    % ws1 = 0.4;
    % Частота среза 2-го ФНЧ
    ws2 = 0.001;
    % ws2 = 0.03;
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

    % yr_double = filter(bandpass_fractional(:,1), 1, adc_input_double(:,1));
    % yr_double = [yr_double(del_proc+1:end); zeros(del_proc,1)]; 
    % 
    % figure(3)
    % plot([yr_double(1:100), adc_input_double(1:100,2)]);
    % 
    % a = adc_input_double(:,1) - yr_double;
    % figure(4)
    % snr(a, sim_options.Fs_sub_adc);


    % snr(double(adc_input_offset(:,1)), sim_options.Fs_sub_adc);

    % запись коэффициентов фильтра дробной задержки в файл 
    % for i = 1:sim_options.M-1
    %     writematrix(coeff_frac_int(:,i), ['src/width_txt/Коэффициенты_фильтров_дробной_задержки_АЦП' num2str(i), '.txt']);
    % end
    %% Синтезируем полосовой фильтр для АЦП0
    lowpass1_adc0 = ws1.*sinc(ws1.*(n'-del_proc));
    lowpass2_adc0 = ws2.*sinc(ws2.*(n'-del_proc));
   
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
    [y3, f3] = freqz(bandpass_fractional(:,4), 1,1024, 'whole', sim_options.Fs_sub_adc);
    for i = 1:sim_options.M
        [y4(:,i), f4(:,i)] = freqz(bandpass_fractional(:,i), 1,1024, 'whole', sim_options.Fs_sub_adc);
    end

    figure(3);
    plot(f4(:,1), abs(y4(:,1)), f4(:,1), abs(y4(:,2)), f4(:,1), abs(y4(:,3)));
    title('АЧХ полосовых фильтров')
    xlabel('Частота') 
    ylabel('Коэффициент передачи') 
    legend({'0.25*Fs','0.5*Fs', '0.75*Fs', 'Эталонный фильтр'},'Location','northeast');

    x2 = xline(500000000, '--', 'Fs/2')
    x2.LabelHorizontalAlignment = 'center'
    x2.LabelVerticalAlignment = 'middle';

    figure(4);
    % fvtool(bandpass_fractional(:,1),bandpass_fractional(:,2));
    plot(f4(:,1), angle(y4(:,1)), f4(:,1), angle(y4(:,2)), f4(:,1), angle(y4(:,3)));
    % fvtool(bandpass_fractional(:,1),bandpass_fractional(:,2));

    %% Расчет коэффициентов фильтра Гилберта

    lowpass_hilbert1 = (2./((n-del_proc)*pi)).*(sin(ws1.*((n-del_proc)*pi)./2)).^2;
    lowpass_hilbert1(1) = 0;
    lowpass_hilbert1(37) = 0;
    lowpass_hilbert2 = (2./((n-del_proc)*pi)).*(sin(ws2.*((n-del_proc)*pi)./2)).^2;
    lowpass_hilbert2(1) = 0;
    lowpass_hilbert2(37) = 0;

    weight_lowpass_hilbert1 = lowpass_hilbert1 .* w_blackman; 
    weight_lowpass_hilbert2 = lowpass_hilbert2 .* w_blackman; 
    bandpass_hilbert = weight_lowpass_hilbert1'  - weight_lowpass_hilbert2';
    % 
    % % Перевод коэффициентов фильтра Гилберта в инты
    hilbert_coeff_int = cast(bandpass_hilbert*2^(sim_options.hilbert_coeff_width-1), sim_options.int_size);
    % % запись коэффициентов фильтра Гилберта в файл
    % writematrix(hilbert_coeff_int, ['src/width_txt/Коэффициенты_фильтра_Гилберта.txt']);

    [y, f] = freqz(double(hilbert_coeff_int)*2^-(sim_options.hilbert_coeff_width-1), 1,1024, 'whole', 1000000000);
    [y1, f1] = freqz(bandpass_hilbert, 1,1024, 'whole', 1000000000);

    % figure(4);
    % plot(f3, abs(y3), f4(:,4), abs(y4(:,4)), f4(:,1), abs(y4(:,1)));
    % title('АЧХ полосовых фильтров дробной задержки и фильтра Гилберта')
    % xlabel('Частота') 
    % ylabel('Коэффициент передачи') 
    % legend('Location','northeast');
   
    figure(4);
    subplot(2,1,1)
    plot(f, abs(y), f, abs(y1));
    title('АЧХ полосовых фильтров дробной задержки и фильтра Гилберта')
    xlabel('Частота') 
    ylabel('Коэффициент передачи') 
    legend('Location','northeast');
    subplot(2,1,2)
    plot(f, angle(y), f, angle(y1));
    %% Коэффициенты эталонного фильтра
    % golden_h = (fir1(128,[0.030 0.440],"bandpass"))';

    %% Данные для переноса сигнала в различные зоны Найквиста
    nn1 = nn' + delay_adc;
    a1 = 2*pi*nn1*Nbp;
    cosi = cos(a1);
    sini = sin(a1);

                                                                %% Первая часть алгоритма калибровки
    %% Эталонный сигнал

    % double
    yr_double = filter(bandpass_adc0, 1, adc_input_double(:,1));
    yr_double = [yr_double(del_proc+1:end); zeros(del_proc,1)]; 

    % adc_input_offset(:,1) = yr_double;
    %% Int
    filter_max_width_out_golden = struct;

    avg = 5*1;

    adc_input_offset(:,1) = double(adc_input(:,1));

    % for k = 1:length(adc_input(:,1))/avg
    %     indexx = avg*(k-1)+1:(avg*k);
    %     y_average_out1 = sum(double(adc_input(indexx,1)));
    %     offset = round(y_average_out1 ./ avg);
    %     adc_input_offset(:,1) = double(adc_input(:,1)) - offset;
    % end
    figure(5);
    subplot(3,1,1)
    plot([adc_input(:,1), adc_input_offset(:,1)]);
    subplot(3,1,2)
    sfdr(double(adc_input(:,1)));
    subplot(3,1,3)
    sfdr(adc_input_offset(:,1));


    % % a = resample(double(adc_input_offset(:,1)),6,5);
    % 
    % interpolate_signal = interp(double(adc_input_offset(:,1)),2);
    % % a =  downsample(double(adc_input_offset(:,1)),2);
    % 
    % figure(6)
    % subplot(3,1,1)
    % plot([double(adc_input_offset(1:100,1)), interpolate_signal(1:100)]);
    % subplot(3,1,2)
    % snr(interpolate_signal, sim_options.Fs_sub_adc);
    % subplot(3,1,3)
    % snr(double(adc_input_offset(:,1)), sim_options.Fs_sub_adc);
    % 
    % for i = 1:3
    %     yri1(:,i) = filter(bandpass_fractional(:,i), 1, interpolate_signal); % filter (стр.6 (15))
    %     % убираем переходной процесс
    %     yri1(:,i) = [yri1(del_proc+1:end,i); zeros(del_proc,1)]; % убираем переходной процесс
    % 
    %     figure(7)
    %     subplot(2,1,1)
    %     plot([double(adc_input(1:100,i)), yri1(1:100,i)]);
    %     subplot(2,1,2)
    %     snr(yri1(:,1), sim_options.Fs_sub_adc*2);
    % 
    %     b(:,i) = downsample(yri1(:,i),2);
    %     figure(8)
    %     subplot(3,1,1)
    %     plot([double(adc_input(1:100,i)), b(1:100,i)]);
    %     subplot(3,1,2)
    %     snr(b(:,i), sim_options.Fs_sub_adc);
    %     subplot(3,1,3)
    %     snr(double(adc_input_offset(:,1)), sim_options.Fs_sub_adc);
    % 
    % end
    % 
    % figure(10)
    % subplot(2,1,1)
    % plot([adc_input(1:200,1), b(1:200,1), b(1:200,2), b(1:200,3)]);
    % subplot(2,1,2)
    % plot([adc_input(1:200,1), adc_input(1:200,2), adc_input(1:200,3), adc_input(1:200,4)]);
    % 
    % % Собираем полученные сигналы в массив для удобства
    % for i = 1:sim_options.M
    %     if i == 1
    %         yri_cut(:,1) = double(adc_input_offset(1:end-del_proc,1));
    %     else
    %         yri_cut(:,i) = b(1:end-del_proc,i-1);
    %     end
    % end
    % sig_adc = zeros(sim_options.M*length(yri_cut(:,1)),1);
    % for i = 1:sim_options.M
    %     sig_adc(i:sim_options.M:end) = yri_cut(:,i);
    % end
    % 
    % figure(9);
    % subplot(3,1,1)
    % snr(s_to_subadc_int, sim_options.Fs/sim_options.Inter);
    % subplot(3,1,2)
    % snr(s_after_subadc, sim_options.Fs/sim_options.Inter);
    % subplot(3,1,3);
    % snr(sig_adc, sim_options.Fs/sim_options.Inter);




	%% Дробная задержка отсчетов сигнала с выхода АЦП0
    for i = 1:sim_options.M-1
        % Фильтр дробной задержки (Int)
        fractional_max_width = readtable(['src/width_txt/Разрядность_максимальных_значений_фильтра_дробной_задержки_АЦП_№' num2str(i) '.xlsx']); 
        %% Фильтр Гилберта (Int)
        hilbert_max_width = readtable(['src/width_txt/Разрядность_максимальных_значений_фильтра_Гилберта_АЦП_№' num2str(i) '.xlsx']);

        [yri, yri_round, ymi, ymi_round, ...
            y_fractional_outInt_div, filter_max_width_out_fractional, ...
            ymi_HilbertInt_div, filter_max_width_out_hilbert] = fractional_part ...
            ( bandpass_fractional(:,i), bandpass_hilbert.', ...
            coeff_frac_int(:,i), hilbert_coeff_int, ...
            adc_input_offset(:,1), ...
            double(adc_input_offset(:,1)), del_proc, fractional_max_width, hilbert_max_width, sim_options);

        % y_fractional_outInt_div = int16(downsample(double(y_fractional_outInt_div),2));
        % ymi_HilbertInt_div = int16(downsample(double(ymi_HilbertInt_div),2));
        % ymi_round = downsample(ymi_round,2);
        % yri_round = downsample(yri_round,2);
        % yri = downsample(yri,2);
        % ymi = downsample(ymi,2);

        %% Ищем макс. значения для записи в файл
        % for j = 1:length(y_fractional_outInt_div)
        %     % возвращаем поделенные положительные значения 
        %     if y_fractional_outInt_div(j) < 0
        %         y_fractional_outInt_abs(j) = y_fractional_outInt_div(j) * cast(-1, sim_options.type_fir_out); % находим число по модулю
        %     else
	    %         y_fractional_outInt_abs(j) = y_fractional_outInt_div(j);
        %     end
        % 
        %     % ищем максимум
        %     if (y_fractional_outInt_abs_max(i) < y_fractional_outInt_abs(j))
        %         y_fractional_outInt_abs_max(i) = y_fractional_outInt_abs(j);
        %     end
        % end

        %% 
        y_fractional_outInt_double = double(y_fractional_outInt_div);
        relative_error_fractional = yri_round./y_fractional_outInt_double;

        figure(6);
        subplot(5,1,1)
        plot([adc_input_offset(1:500,1), yri(1:500)]);
        title('Выходные сигналы фильтров дробной задержки');
        xlabel('Номер отсчета');
        ylabel('Значение отсчета'); 
        legend({'double','int'},'Location','northeast');
        subplot(5,1,2);
        snr(yri, sim_options.Fs_sub_adc/100);
        subplot(5,1,3);
        snr(yri_round, sim_options.Fs_sub_adc);
        subplot(5,1,4);
        snr(y_fractional_outInt_double, sim_options.Fs_sub_adc);
        subplot(5,1,5);
        plot(relative_error_fractional);
        title('Относительная ошибка выходного сигнала фильтра дробной задержки между double и integer');
        xlabel('Номер отсчета');
        ylabel('Значение ошибки');
        x4 = xline(37, '--', 'Переходной процесс фильтра');
        x4.LabelHorizontalAlignment = 'center';
        x4.LabelVerticalAlignment = 'middle';

        % %% Ищем макс. значения для записи в файл
        % for j = 1:length(ymi_HilbertInt_div)
        %     % возвращаем поделенные положительные значения 
        %     if ymi_HilbertInt_div(j) < 0
        %         ymi_HilbertInt_abs(j) = ymi_HilbertInt_div(j) * cast(-1, sim_options.type_fir_out); % находим число по модулю
        %     else
	    %         ymi_HilbertInt_abs(j) = ymi_HilbertInt_div(j);
        %     end
        % 
        %     % ищем максимум
        %     if (ymi_HilbertInt_abs_max(i) < ymi_HilbertInt_abs(j))
        %         ymi_HilbertInt_abs_max(i) = ymi_HilbertInt_abs(j);
        %     end
        % end

        %%
        ymi_HilbertInt_double = double(ymi_HilbertInt_div);
        relative_error_hilbert = ymi_round./ymi_HilbertInt_double;

        figure(7);
        subplot(5,1,1)
        plot([ymi_round(1:500), double(ymi_HilbertInt_div(1:500))]);
        title('Выходные сигналы фильтра Гилберта');
        xlabel('Номер отсчета'); 
        ylabel('Значение отсчета');
        legend({'double','int'},'Location','northeast');
        subplot(5,1,2)
        snr(ymi, 1000000000);
        subplot(5,1,3);
        snr(ymi_round, 1000000000);
        subplot(5,1,4)
        snr(double(ymi_HilbertInt_div), 1000000000);
        subplot(5,1,5)
        plot(relative_error_hilbert);
        title('Относительная ошибка выходного сигнала фильтра Гилберта между double и integer');
        xlabel('Номер отсчета'); 
        ylabel('Значение ошибки'); 
        x4 = xline(73, '--', 'Переходной процесс фильтра');
        x4.LabelHorizontalAlignment = 'center';
        x4.LabelVerticalAlignment = 'middle';

        %% Перенос сигнала в заданные зоны Найквиста
        [yric(:,i), yric_int(:,i)] = single_sideband(yri, y_fractional_outInt_div, ymi, ymi_HilbertInt_div, cosi(:,i), sini(:,i), i, sim_options);

    end

    yri_cut = zeros(length(adc_input_offset(1:end-del_proc,1)),sim_options.M);
    yri_cut_int = cast(zeros(length(adc_input_offset(1:end-del_proc,1)),sim_options.M), sim_options.type_fir_out);

    % Собираем полученные сигналы в массив для удобства
    for i = 1:sim_options.M
        if i == 1
            yri_cut(:,1) = double(adc_input_offset(1:end-del_proc,1));
            yri_cut_int(:,1) = adc_input_offset(1:end-del_proc,1);
        else
            yri_cut(:,i) = yric(1:end-del_proc,i-1);
            yri_cut_int(:,i) = yric_int(1:end-del_proc,i-1);
            % yri_cut(:,i) = double(adc_input(1:2994-del_proc,i));
            % yri_cut_int(:,i) = double(adc_input(1:2994-del_proc,i));
        end
    end

    % figure(10)
    % subplot(2,1,1)
    % plot([yri_cut(1:200,1), yri_cut(1:200,2), yri_cut(1:200,3), yri_cut(1:200,4)]);
    % subplot(2,1,2)
    % plot([adc_input(1:200,1), adc_input(1:200,2), adc_input(1:200,3), adc_input(1:200,4)]);

    %% Тестируем первую часть алгоритма калибровки

    % sig_adc_gold = zeros(sim_options.M*length(adc_input_offset(:,1)),1);
    sig_adc = zeros(sim_options.M*length(yri_cut(:,1)),1);
    sig_adc_int = zeros(sim_options.M*length(yri_cut_int(:,1)),1);
	for i = 1:sim_options.M
        % sig_adc_gold(i:sim_options.M:end) = adc_input_offset(:,i);
        sig_adc(i:sim_options.M:end) = yri_cut(:,i);
        sig_adc_int(i:sim_options.M:end) = double(yri_cut_int(:,i));
    end

    %%
    % save (sprintf(num2str(clock) + ".mat"));
    % load ('2025             12             29             10              5         25.901.mat'); 
    % load ('2026              1             23             16             12         35.392.mat'); 
    % sim_options.Ls = 500;
    
    % load ('test_after_fractional_filters_50_MHz'); % 317
    % sim_options.enable_log = false;
    % sim_options.type_2x2_det = "int64";
    % sim_options.type_3x3_det = "double";
    % sim_options.type_4x4_det = "double";
    % sim_options.type_5x5_det = "double";
    %%
    figure(8);
    subplot(3,1,1)
    plot([s_to_subadc_int(1:750), sig_adc_int(1:750)]);
    title('Исходный сигнал');
    xlabel('Номер отсчета'); 
    ylabel('Амплитуда')
    subplot(3,1,2)
    plot(sig_adc(1:750));
    title('Сигнал после фильтров double и int');
    xlabel('Номер отсчета'); 
    ylabel('Амплитуда')
    subplot(3,1,3)
    plot(sig_adc ./ sig_adc_int);
    title('Относительная ошибка double и int');
    xlabel('Номер отсчета'); 
    ylabel('Величина ошибки')

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

        % for k = 1:length(adc_input(:,1))/avg
        %     indexx = avg*(k-1)+1:(avg*k);
        %     y_average_out = sum(double(adc_input(indexx,i+1)));
        %     offset = y_average_out ./ avg;
        %     adc_input_offset(:,i+1) = double(adc_input(:,i+1)) - offset;
        % end
        % 
        % figure(3)
        % subplot(5,1,1)
        % plot([double(adc_input(1:100,i+1)), adc_input_offset(1:100,i+1)]);
        % subplot(5,1,2)
        % snr(double(adc_input(:,i+1)));
        % subplot(5,1,3)
        % sfdr(double(adc_input(:,i+1)));
        % subplot(5,1,4)
        % snr(adc_input_offset(:,i+1));
        % subplot(5,1,5)
        % sfdr(adc_input_offset(:,i+1));

        adc_input_offset(:,i+1) = double(adc_input(:,i+1));

        % det_max_width = readtable(['src/width_txt/Разрядность_максимальных_значений_определителя_АЦП_№' num2str(i) '.xlsx']); 
        % read_max_width.det_max_width = table2array(det_max_width(:,2:end));
        % 
        % adaptive_max_width = readtable(['src/width_txt/Разрядность_максимальных_значений_адаптивного_фильтра_АЦП_№' num2str(i) '.xlsx']); 
        % read_max_width.adaptive_max_width = (adaptive_max_width(:,2:end));


        % [y_average_outInt(:,i+1), filter_max_width_out_average ...
        % ] = fir_filter(coeff_average_int, adc_input(:,i+1), 2, golden_max_width, sim_options.width_fractional, sim_options);
        % y_average_outInt_del(:,i+1) = int16(bitshift(y_average_outInt(:,i+1),-2));


        [y_array(:,i), y_array_double(:,i), y_array_int(:,i), determinate_struct(i), Divide_max(i), adaptive_filter_struct_max(i) ...
            ] = least_mean_squares(adc_input_double(:,i+1), adc_input_offset(:,i+1), yri_cut(:,i+1), yri_cut_int(:,i+1), read_max_width, sim_options);
        
   end

    %% create main signal after LS algorithm (switch after sub-adc)
    x_after_adc = zeros(length(y_array)*sim_options.M,1);
    x_after_adc_double = zeros(length(y_array_double)*sim_options.M,1);
    x_after_adc_int = zeros(length(y_array_int)*sim_options.M,1);
    
    for i = 1:sim_options.M
        if i == 1
            x_after_adc(i:sim_options.M:end) = adc_input_offset(1:length(y_array),1);
            x_after_adc_double(i:sim_options.M:end) = double(yri_cut_int(1:length(y_array_double),1));
            x_after_adc_int(i:sim_options.M:end) = double(yri_cut_int(1:length(y_array_int),1));
        else
            x_after_adc(i:sim_options.M:end) = y_array(:,i-1);
            x_after_adc_double(i:sim_options.M:end) = y_array_double(:,i-1); % * 2^-(sim_options.divide_factor);
            x_after_adc_int(i:sim_options.M:end) = double(y_array_int(:,i-1)); % * 2^-(sim_options.divide_factor);
        end
    end

    % figure(20);
    % plot(double(sig_adc_int(1:3500) ./ x_after_adc_int(1:3500)))
    % title('Сходимость алгоритма калибровки');
    % xlabel('Количество отсчетов'); 
    % ylabel('Результат относительной ошибки'); 
    % 
    figure(11);
    subplot(5,1,1);
    plot([s_to_subadc_int(1:600), s_after_subadc(1:600), x_after_adc(1:600)]);
    title('Исходный сигнал и выход алгоритма LU');
    xlabel('Номер отсчета'); 
    ylabel('Амплитуда'); 
    legend({'Исходный сигнал','Сигнал с выхода адаптивного фильтра double'},'Location','northeast');
    subplot(5,1,2);
    plot([s_to_subadc_int(1:1000), x_after_adc(1:1000)]);
    title('Исходный сигнал и выход адаптивного фильтра matlab');
    xlabel('Номер отсчета'); 
    ylabel('Амплитуда'); 
    legend({'Исходный сигнал','Сигнал с выхода адаптивного фильтра '},'Location','northeast');
    subplot(5,1,3);
    plot([s_to_subadc_int(1:1000), x_after_adc_double(1:1000)]);
    title('Исходный сигнал и выход адаптивного фильтра double');
    xlabel('Номер отсчета'); 
    ylabel('Амплитуда'); 
    legend({'Исходный сигнал','Сигнал с выхода адаптивного фильтра '},'Location','northeast');
    subplot(5,1,4);
    plot([s_to_subadc_int(1:1000), x_after_adc_int(1:1000)]);
    title('Исходный сигнал и выход адаптивного фильтра int');
    xlabel('Номер отсчета'); 
    ylabel('Амплитуда');
    legend({'Исходный сигнал','Сигнал с выхода адаптивного фильтра int'},'Location','northeast');
    subplot(5,1,5);
    plot([double(s_to_subadc_int(1:length(x_after_adc_double))) ./ x_after_adc_double]);
    %% SFDR
    figure(12);
    subplot(5,1,1);
    sfdr(s_to_subadc_int(1:length(x_after_adc)), sim_options.Fs/sim_options.Inter);
    subplot(5,1,2);
    sfdr(s_after_subadc(1:length(x_after_adc)), sim_options.Fs/sim_options.Inter);
    subplot(5,1,3);
    sfdr(x_after_adc(1:length(x_after_adc)), sim_options.Fs/sim_options.Inter);
    subplot(5,1,4);
    sfdr(x_after_adc_double(1:length(x_after_adc_double)), sim_options.Fs/sim_options.Inter);
    subplot(5,1,5);
    sfdr(x_after_adc_int(1:length(x_after_adc)), sim_options.Fs/sim_options.Inter);
    %% SNR
    figure(13);
    subplot(5,1,1);
    snr(s_to_subadc_int(1:length(x_after_adc)), sim_options.Fs/sim_options.Inter);
    subplot(5,1,2);
    snr(s_after_subadc(1:length(x_after_adc)), sim_options.Fs/sim_options.Inter);
    subplot(5,1,3);
    snr(x_after_adc(1:length(x_after_adc)), sim_options.Fs/sim_options.Inter);
    subplot(5,1,4);
    snr(x_after_adc_double(1:length(x_after_adc_double)), sim_options.Fs/sim_options.Inter);
    subplot(5,1,5);
    snr(x_after_adc_int(1:length(x_after_adc)), sim_options.Fs/sim_options.Inter);

end
