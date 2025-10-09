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

freq = sim_options.freq;

for num = 1:sim_options.num_cycles
 
    Z = ceil(freq/(sim_options.Fs/sim_options.Inter/2/sim_options.M));      % Nyquist zone

    [s_to_subadc, s_to_subadc_int, adc_input, adc_input_int, s_after_subadc, s_after_subadc_int] = gen_oversampled_signal(sim_options.M, sim_options.Fs, freq, sim_options.SNR, ...
        sim_options.Inter, sim_options.StopTime, sim_options.MODEL_ERROR, sim_options.time_skew_array, sim_options.gain_error_array);

    [x_after_adc, x_after_adc_int, snr_s] = adc_calibration(sim_options, adc_input_int, s_to_subadc_int, s_after_subadc, Z);

    %% Measurements1
    figure(5);
    subplot(2,1,1)
    plot([x_after_adc(1:500)])
    title('Исходный сигнал до искажения и выход адаптивного фильтра (double)')
    xlabel('Номер отсчета') 
    ylabel('Амплитуда') 
    legend({'Исходный сигнал','Сигнал с выхода адаптивного фильтра'},'Location','northeast')

    subplot(2,1,2)
    plot([s_to_subadc_int(1:length(x_after_adc_int)), x_after_adc_int]); %, error_out(:,2), error_out(:,3)]);
    title('Исходный сигнал до искажения и выход адаптивного фильтра (int)')
    xlabel('Номер отсчета') 
    ylabel('Отношение') 
    legend({'Исходный сигнал','Сигнал с выхода адаптивного фильтра'},'Location','northeast')
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
    % 
    figure(7);
    subplot(4,1,1);
    snr(s_to_subadc(1:length(x_after_adc)), sim_options.Fs/sim_options.Inter);
    subplot(4,1,2);
    snr(s_after_subadc(1:length(x_after_adc)), sim_options.Fs/sim_options.Inter);
    subplot(4,1,3);
    snr(x_after_adc(1:length(x_after_adc)), sim_options.Fs/sim_options.Inter);
    subplot(4,1,4);
    snr(x_after_adc_int(1:length(x_after_adc)), sim_options.Fs/sim_options.Inter);
    % 
    % snr_in_double(num) = snr(s_after_subadc, sim_options.Fs/sim_options.Inter);
    % snr_in_int(num) = snr(double(s_after_subadc_int), sim_options.Fs/sim_options.Inter);
    % snr_output_double(num) = snr(x_after_adc, sim_options.Fs/sim_options.Inter);
    % snr_output_int(num) = snr(x_after_adc_int, sim_options.Fs/sim_options.Inter);
    % 
    % sfdr_in_double(num) = sfdr(s_after_subadc, sim_options.Fs/sim_options.Inter);
    % sfdr_in_int(num) = sfdr(double(s_after_subadc_int), sim_options.Fs/sim_options.Inter);
    % sfdr_output_double(num) = sfdr(x_after_adc, sim_options.Fs/sim_options.Inter);
    % sfdr_output_int(num) = sfdr(x_after_adc_int, sim_options.Fs/sim_options.Inter);
    
    norm_freq(num) = freq/(sim_options.Fs/sim_options.Inter/sim_options.M);
    num_array(:,num) = num;


    % freq
    freq = freq + sim_options.step; % frequency of fundamental tone

    % SNR
    sim_options.SNR = sim_options.SNR + sim_options.Step_of_SNR;
    snr_double(num) = snr_s(1);
    snr_16(num) = snr_s(2);
    snr_19(num) = snr_s(3);
    snr_21(num) = snr_s(4);
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

    figure(9);
    plot(norm_freq, snr_double, '-o', norm_freq, snr_16, '-o', norm_freq, snr_19, '-o', norm_freq, snr_21, '-o');
    title('Зависимость разрядности коэффициентов на выходной итоговый сигнал')
    xlabel('Нормированная частота') 
    ylabel('SNR (dB)') 
    legend({'double', '16 бит', '19 бит', '21 бит'}, 'Location','northwest');


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