function runsim(sim_options)

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

for num = 1:sim_options.num_cycles

    sim_options.SNR = sim_options.SNR + 5;
    snr_array(num) = sim_options.SNR; 

    [s_to_subadc, adc_input, adc_input_int, s_after_subadc, sim_options] = gen_oversampled_signal(sim_options);
    [sig_adc, x_after_adc, x_after_adc_int, error_det, error_det_lu] = adc_calibration(sim_options, adc_input, adc_input_int, s_to_subadc, s_after_subadc);

    %% Measurements1
    figure(4);
    % subplot(2,1,1)
    plot([s_to_subadc(1:length(x_after_adc)), x_after_adc])
    % title('Отношение между отсчетами I-составляющей')
    xlabel('Номер отсчета') 
    ylabel('Амплитуда') 
    legend({'Исходный сигнал','Сигнал с выхода адаптивного фильтра'},'Location','northeast')
    % subplot(2,1,2)
    % plot([error_out(:,1)]); %, error_out(:,2), error_out(:,3)]);
    % title('Относительная ошибка между исходным сигналом и выходом адаптивного фильтра')
    % xlabel('Номер отсчета') 
    % ylabel('Отношение') 

    figure(5);
    subplot(4,1,1);
    sfdr(s_to_subadc(1:length(x_after_adc)), sim_options.Fs/sim_options.Inter);
    subplot(4,1,2);
    sfdr(s_after_subadc(1:length(x_after_adc)), sim_options.Fs/sim_options.Inter);
    subplot(4,1,3);
    sfdr(x_after_adc(1:length(x_after_adc)), sim_options.Fs/sim_options.Inter);
    subplot(4,1,4);
    sfdr(x_after_adc_int(1:length(x_after_adc)), sim_options.Fs/sim_options.Inter);

    figure(6);
    subplot(4,1,1);
    snr(s_to_subadc(1:length(x_after_adc)), sim_options.Fs/sim_options.Inter);
    subplot(4,1,2);
    snr(s_after_subadc(1:length(x_after_adc)), sim_options.Fs/sim_options.Inter);
    subplot(4,1,3);
    snr(x_after_adc(1:length(x_after_adc)), sim_options.Fs/sim_options.Inter);
    subplot(4,1,4);
    snr(x_after_adc_int(1:length(x_after_adc)), sim_options.Fs/sim_options.Inter);

    snr_in_id(num) = snr(s_after_subadc, sim_options.Fs/sim_options.Inter);
    snr_input(num) = snr(x_after_adc_int, sim_options.Fs/sim_options.Inter);
    snr_output(num) = snr(x_after_adc, sim_options.Fs/sim_options.Inter);

    sfdr_in_id(num) = sfdr(s_after_subadc, sim_options.Fs/sim_options.Inter);
    sfdr_input(num) = sfdr(x_after_adc_int, sim_options.Fs/sim_options.Inter);
    sfdr_output(num) = sfdr(x_after_adc, sim_options.Fs/sim_options.Inter);
    norm_freq(num) = sim_options.freq/(sim_options.Fs/sim_options.Inter/sim_options.M);


    error_det_array(:,num) = error_det; 
    error_det_array_lu(:,num) = error_det_lu;
    num_array(:,num) = num;

end

    figure(7);
    subplot(2,1,1)
    plot(num_array, snr_input, '-o', num_array, snr_output, '-o', num_array, snr_in_id, '-o');
    title('SNR')
    xlabel('Номер итерации') 
    ylabel('SNR (dB)') 
    legend('Определитель 70 бит', 'double', 'Исходный сигнал с ошибками')
    subplot(2,1,2)
    plot(num_array, sfdr_input, '-o', num_array, sfdr_output, '-o', num_array, sfdr_in_id, '-o');
    title('SFDR (dB)')
    xlabel('Номер итерации') 
    ylabel('SFDR (dB)') 
    legend('Определитель 70 бит', 'double', 'Исходный сигнал с ошибками')


    % figure(8);
    % % subplot(7,1,1)
    % plot([error_det_array(:,1), error_det_array_lu(:,1)]);
    % title('Относительная ошибка определителей при SNR 10')
    % xlabel('Номер отсчета') 
    % ylabel('Величина ошибки') 
    % subplot(7,1,2)
    % plot([error_det_array(:,2), error_det_array_lu(:,2)]);
    % title('Относительная ошибка определителей при SNR 20')
    % xlabel('Номер отсчета') 
    % ylabel('Величина ошибки')
    % subplot(7,1,3)
    % plot([error_det_array(:,3), error_det_array_lu(:,3)]);
    % title('Относительная ошибка определителей при SNR 30')
    % xlabel('Номер отсчета') 
    % ylabel('Величина ошибки')
    % subplot(7,1,4)
    % plot([error_det_array(:,4), error_det_array_lu(:,4)]);
    % title('Относительная ошибка определителей при SNR 40')
    % xlabel('Номер отсчета') 
    % ylabel('Величина ошибки')
    % subplot(7,1,5)
    % plot([error_det_array(:,5), error_det_array_lu(:,5)]);
    % title('Относительная ошибка определителей при SNR 50')
    % xlabel('Номер отсчета') 
    % ylabel('Величина ошибки')
    % subplot(7,1,6)
    % plot([error_det_array(:,6), error_det_array_lu(:,6)]);
    % title('Относительная ошибка определителей при SNR 60')
    % xlabel('Номер отсчета') 
    % ylabel('Величина ошибки')
    % subplot(7,1,7)
    % plot([error_det_array(:,7), error_det_array_lu(:,7)]);
    % title('Относительная ошибка определителей при SNR 70')
    % xlabel('Номер отсчета') 
    % ylabel('Величина ошибки')

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