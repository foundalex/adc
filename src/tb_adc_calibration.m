function tb_adc_calibration (sim_options)

all_figs = findobj(0, 'type', 'figure');
delete(setdiff(all_figs, 1));
clc;

% Set Random number generators initial state
% reset random number generators based on current clock value
rand('state',sum(100*clock));
randn('state',sum(100*clock));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main simulation loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize simulation timer
start_time = clock;

%%
width_mult = int8(zeros(sim_options.N,sim_options.M-1));
width_sum = int8(zeros(sim_options.N-1,sim_options.M-1));

fractional_mult_max = cast(zeros(sim_options.N,sim_options.M-1), sim_options.int_size);
fractional_sum_max = cast(zeros(sim_options.N-1,sim_options.M-1), sim_options.int_size);

fractional_total_width_mult_max = cast(zeros(sim_options.N,sim_options.M-1), sim_options.int_size);
fractional_total_width_sum_max = cast(zeros(sim_options.N-1,sim_options.M-1), sim_options.int_size);

%%
width_mult_h = int8(zeros(sim_options.N, sim_options.M-1));
width_sum_h = int8(zeros(sim_options.N-1, sim_options.M-1));

hilbert_mult_max = cast(zeros(sim_options.N,sim_options.M-1), sim_options.int_size);
hilbert_sum_max = cast(zeros(sim_options.N-1,sim_options.M-1), sim_options.int_size);

hilbert_total_width_mult_max = cast(zeros(sim_options.N,sim_options.M-1), sim_options.int_size);
hilbert_total_width_sum_max = cast(zeros(sim_options.N-1,sim_options.M-1), sim_options.int_size);

%%
for num = 1:sim_options.num_cycles
 
    Z = ceil(sim_options.freq/(sim_options.Fs/sim_options.Inter/2/sim_options.M));      % Nyquist zone

    [s_to_subadc, s_to_subadc_int, adc_input, adc_input_int, s_after_subadc, s_after_subadc_int] = gen_oversampled_signal(sim_options.M, sim_options.Fs, ...
        sim_options.freq, sim_options.SNR, sim_options.Inter, sim_options.StopTime, sim_options.MODEL_ERROR, sim_options.time_skew_array, sim_options.gain_error_array);

    [x_after_adc, x_after_adc_int, snr_s, fractional_mult, fractional_sum, fractional_width_total_mult, fractional_width_total_sum, ...
        hilbert_width_mult, hilbert_width_sum, hilbert_width_total_mult, hilbert_width_total_sum, ...
        ... % Determinant
        Det2x2_mult1_abs, Det2x2_mult2_abs, Det2x2_sum_abs] ...
        = adc_calibration(sim_options, adc_input_int, s_to_subadc_int, s_after_subadc, Z);

    % Записываем значения каждого фильтра
    for i = 1:sim_options.M-1
        for j = 1:sim_options.N
            % выбираем максимальное значение сигнала умножителей фильтра
            % дробной задержки
            if (fractional_mult_max(j,i) < fractional_mult(j,i))
                fractional_mult_max(j,i) = fractional_mult(j,i); 
            end

            % выбираем максимальное значение разрядности умножителя фильтра
            % дробной задержки
            if (fractional_total_width_mult_max(j,i) < fractional_width_total_mult(j,i))
                fractional_total_width_mult_max(j,i) = fractional_width_total_mult(j,i); 
            end
            %%
            % выбираем максимальное значение сигнала умножителей фильтра
            % Гилберта
            if (hilbert_mult_max(j,i) < hilbert_width_mult(j,i))
                hilbert_mult_max(j,i) = hilbert_width_mult(j,i); 
            end

            % выбираем максимальное значение разрядности умножителя фильтра
            % Гилберта
            if (hilbert_total_width_mult_max(j,i) < hilbert_width_total_mult(j,i))
                hilbert_total_width_mult_max(j,i) = hilbert_width_total_mult(j,i); 
            end

        end


        %%
        for j = 1:sim_options.N-1
            % выбираем максимальное значение сигнала сумматоров
            % фильтра дробной задержки 
            if (fractional_sum_max(j,i) < fractional_sum(j,i))
                fractional_sum_max(j,i) = fractional_sum(j,i); 
            end

            % выбираем максимальное значение разрядности сумматоров
            % фильтра дробной задержки
            if (fractional_total_width_sum_max(j,i) < fractional_width_total_sum(j,i))
                fractional_total_width_sum_max(j,i) = fractional_width_total_sum(j,i); 
            end

            %%
            % выбираем максимальное значение сигнала сумматоров
            % фильтра Гилберта
            if (hilbert_sum_max(j,i) < hilbert_width_sum(j,i))
                hilbert_sum_max(j,i) = hilbert_width_sum(j,i); 
            end

            % выбираем максимальное значение разрядности сумматоров
            % фильтра Гилберта
            if (hilbert_total_width_sum_max(j,i) < hilbert_width_total_sum(j,i))
                hilbert_total_width_sum_max(j,i) = hilbert_width_total_sum(j,i); 
            end

        end

    end

    %% Measurements1
    % figure(5);
    % subplot(2,1,1)
    % plot([x_after_adc(1:500)])
    % title('Исходный сигнал до искажения и выход адаптивного фильтра (double)')
    % xlabel('Номер отсчета') 
    % ylabel('Амплитуда') 
    % legend({'Исходный сигнал','Сигнал с выхода адаптивного фильтра'},'Location','northeast')
    % 
    % subplot(2,1,2)
    % plot([s_to_subadc_int(1:length(x_after_adc_int)), x_after_adc_int]); %, error_out(:,2), error_out(:,3)]);
    % title('Исходный сигнал до искажения и выход адаптивного фильтра (int)')
    % xlabel('Номер отсчета') 
    % ylabel('Отношение') 
    % legend({'Исходный сигнал','Сигнал с выхода адаптивного фильтра'},'Location','northeast')
    %%
    % figure(6);
    % subplot(4,1,1);
    % sfdr(s_to_subadc(1:length(x_after_adc)), sim_options.Fs/sim_options.Inter);
    % subplot(4,1,2);
    % sfdr(s_after_subadc(1:length(x_after_adc)), sim_options.Fs/sim_options.Inter);
    % subplot(4,1,3);
    % sfdr(x_after_adc(1:length(x_after_adc)), sim_options.Fs/sim_options.Inter);
    % subplot(4,1,4);
    % sfdr(x_after_adc_int(1:length(x_after_adc)), sim_options.Fs/sim_options.Inter);
    % % 
    % figure(7);
    % subplot(4,1,1);
    % snr(s_to_subadc(1:length(x_after_adc)), sim_options.Fs/sim_options.Inter);
    % subplot(4,1,2);
    % snr(s_after_subadc(1:length(x_after_adc)), sim_options.Fs/sim_options.Inter);
    % subplot(4,1,3);
    % snr(x_after_adc(1:length(x_after_adc)), sim_options.Fs/sim_options.Inter);
    % subplot(4,1,4);
    % snr(x_after_adc_int(1:length(x_after_adc)), sim_options.Fs/sim_options.Inter);
    % % 
    % snr_in_double(num) = snr(s_after_subadc, sim_options.Fs/sim_options.Inter);
    % snr_in_int(num) = snr(double(s_after_subadc_int), sim_options.Fs/sim_options.Inter);
    % snr_output_double(num) = snr(x_after_adc, sim_options.Fs/sim_options.Inter);
    % snr_output_int(num) = snr(x_after_adc_int, sim_options.Fs/sim_options.Inter);
    % 
    % sfdr_in_double(num) = sfdr(s_after_subadc, sim_options.Fs/sim_options.Inter);
    % sfdr_in_int(num) = sfdr(double(s_after_subadc_int), sim_options.Fs/sim_options.Inter);
    % sfdr_output_double(num) = sfdr(x_after_adc, sim_options.Fs/sim_options.Inter);
    % sfdr_output_int(num) = sfdr(x_after_adc_int, sim_options.Fs/sim_options.Inter);
    % 
    % norm_freq(num) = freq/(sim_options.Fs/sim_options.Inter/sim_options.M);
    % num_array(:,num) = num;

    % freq
    sim_options.freq = sim_options.freq + sim_options.step; % frequency of fundamental tone
    % SNR
    sim_options.SNR = sim_options.SNR + sim_options.Step_of_SNR;
end


if sim_options.enable_mask == false
    for i = 1:sim_options.M-1
        for j = 1:length(fractional_mult_max(:,i))
            width_mult(j,i) = define_of_width_int(fractional_mult_max(j,i), sim_options.int_size, sim_options.width_fractional);
            width_mult_h(j,i) = define_of_width_int(hilbert_mult_max(j,i), sim_options.int_size, sim_options.width_hilbert);
        end

        for j = 1:length(fractional_sum_max(:,i))
            width_sum(j,i) = define_of_width_int(fractional_sum_max(j,i), sim_options.int_size, sim_options.width_fractional);
            width_sum_h(j,i) = define_of_width_int(hilbert_sum_max(j,i), sim_options.int_size, sim_options.width_hilbert);
        end
        %% Запись данных для фильтра дробной задержки
        % запись макс. значений сигнала
        writematrix(fractional_mult_max(:,i), ['src/width_txt/Max_value_multiplier_Fractional_filter_' num2str(i) '.txt']);
        writematrix(fractional_sum_max(:,i), ['src/width_txt/Max_value_adder_Fractional_filter_' num2str(i) '.txt']);
        % запись разрядности макс. значений сигнала
        writematrix(width_mult(:,i), ['src/width_txt/Width_multiplier_Fractional_filter_' num2str(i) '.txt']);
        writematrix(width_sum(:,i), ['src/width_txt/Width_adder_Fractional_filter_' num2str(i) '.txt']);
        % запись суммарной разрядности сумматоров и умножителей
        writematrix(fractional_total_width_mult_max(:,i), ['src/width_txt/Total_width_multiplier_Fractional_filter_' num2str(i) '.txt']);
        writematrix(fractional_total_width_sum_max(:,i), ['src/width_txt/Total_width_adder_Fractional_filter_' num2str(i) '.txt']);
        %% Запись данных для фильтра Гилберта
        % запись макс. значений сигнала
        writematrix(hilbert_mult_max(:,i), ['src/width_txt/Max_value_multiplier_Hilbert_filter_' num2str(i) '.txt']);
        writematrix(hilbert_sum_max(:,i), ['src/width_txt/Max_value_adder_Hilbert_filter_' num2str(i) '.txt']);
        % запись разрядности макс. значений сигнала
        writematrix(width_mult_h(:,i), ['src/width_txt/Width_multiplier_Hilbert_filter_' num2str(i) '.txt']);
        writematrix(width_sum_h(:,i), ['src/width_txt/Width_adder_Hilbert_filter_' num2str(i) '.txt']);
        % запись суммарной разрядности сумматоров и умножителей
        writematrix(hilbert_total_width_mult_max(:,i), ['src/width_txt/Total_width_multiplier_Hilbert_filter_' num2str(i) '.txt']);
        writematrix(hilbert_total_width_sum_max(:,i), ['src/width_txt/Total_width_adder_Hilbert_filter_' num2str(i) '.txt']);
    end
end

    % figure(8);
    % subplot(2,1,1)
    % plot(norm_freq, snr_in_double, '-o', norm_freq, snr_in_int, '-o', norm_freq, snr_output_double, '-o', norm_freq, snr_output_int, '-o');
    % title('SNR')
    % xlabel('Нормированная частота') 
    % ylabel('SNR (dB)') 
    % legend({'Входной сигнал c ошибками double', 'Входной сигнал с ошибками int', 'Выходной сигнал double', 'Выходной сигнал int'}, 'Location','northwest');
    % 
    % subplot(2,1,2)
    % plot(norm_freq, sfdr_in_double, '-o', norm_freq, sfdr_in_int, '-o', norm_freq, sfdr_output_double, '-o', norm_freq, sfdr_output_int, '-o');
    % title('SFDR (dB)')
    % xlabel('Нормированная частота') 
    % ylabel('SFDR (dB)') 
    % legend({'Входной сигнал c ошибками double', 'Входной сигнал с ошибками int', 'Выходной сигнал double', 'Выходной сигнал int'}, 'Location','northwest');

    % figure(9);
    % plot(norm_freq, snr_double, '-o', norm_freq, snr_16, '-o', norm_freq, snr_19, '-o', norm_freq, snr_21, '-o');
    % title('Зависимость разрядности коэффициентов на выходной итоговый сигнал')
    % xlabel('Нормированная частота') 
    % ylabel('SNR (dB)') 
    % legend({'double', '16 бит', '19 бит', '21 бит'}, 'Location','northwest');


    %% Measurements2
    % figure(7);
    % subplot(2,1,1)
    % plot(norm_freq, snr_in_id, '-o', norm_freq, snr_input, '-o', norm_freq, snr_output, '-o');
    % title('SNR')
    % xlabel('Нормированная частота') 
    % ylabel('SNR (dB)') 
    % legend('до калибровки без искажений', 'до калибровки с искажениями', 'после калибровки')
    % subplot(2,1,2)
    % plot(norm_freq, sfdr_in_id, '-o', norm_freq, sfdr_input, '-o', norm_freq, sfdr_output, '-o');
    % title('SFDR (dB)')
    % xlabel({'Нормированная частота fнорм = f/(Fs/M)','Fs - частота дискретизации всего TI-ADC, М - количество каналов'}) 
    % ylabel('SFDR (dB)') 
    % legend('до калибровки без искажений', 'до калибровки с искажениями','после калибровки')

    % x4 = xline(0.42, '--', 'Интервал из статьи 1-ой зоны Найквиста')
    % x4.LabelHorizontalAlignment = 'center'
    % x4.LabelVerticalAlignment = 'middle';
    % x2 = xline(0.55, '--', 'Интервал из статьи начало 2-ой зоны Найквиста')
    % x2.LabelHorizontalAlignment = 'center'
    % x2.LabelVerticalAlignment = 'middle';
    % x3 = xline(0.92, '--', 'Интервал из статьи конец 2-ой зоны Найквиста')
    % x3.LabelHorizontalAlignment = 'center'
    % x3.LabelVerticalAlignment = 'middle';
    % y2 = yline(79,'--', 'Нижняя граница SFDR (dB)')
    % y2.LabelHorizontalAlignment = 'left'


stop_time = clock;
elapsed_time = etime(stop_time,start_time);

fprintf('Simulation duration: %g seconds\n',elapsed_time);