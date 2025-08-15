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

    [x_to_subadc, adc_input, x_after_subadc, sim_options] = gen_oversampled_signal(sim_options);
    [sig_adc, x_after_adc, error_out] = adc_calibration(sim_options, adc_input);

    %% Measurements1
    figure(4);
    subplot(2,1,1)
    plot([x_after_subadc(1:150), x_after_adc(1:150)])
    % title('Отношение между отсчетами I-составляющей')
    xlabel('Номер отсчета') 
    ylabel('Амплитуда') 
    legend({'Исходный сигнал','Сигнал с выхода адаптивного фильтра'},'Location','northeast')
    subplot(2,1,2)
    plot([error_out(:,1)]); %, error_out(:,2), error_out(:,3)]);
    title('Относительная ошибка между исходным сигналом и выходом адаптивного фильтра')
    xlabel('Номер отсчета') 
    ylabel('Отношение') 

    figure(5);
    subplot(4,1,1);
    sfdr(x_to_subadc(1:length(x_after_adc)), sim_options.Fs/sim_options.Inter);
    subplot(4,1,2);
    sfdr(x_after_subadc(1:length(x_after_adc)), sim_options.Fs/sim_options.Inter);
    subplot(4,1,3);
    sfdr(sig_adc(1:length(x_after_adc)), sim_options.Fs/sim_options.Inter);
    subplot(4,1,4);
    sfdr(x_after_adc(1:length(x_after_adc)), sim_options.Fs/sim_options.Inter);

    figure(6);
    subplot(4,1,1);
    snr(x_to_subadc(1:length(x_after_adc)), sim_options.Fs/sim_options.Inter);
    subplot(4,1,2);
    snr(x_after_subadc(1:length(x_after_adc)), sim_options.Fs/sim_options.Inter);
    subplot(4,1,3);
    snr(sig_adc(1:length(x_after_adc)), sim_options.Fs/sim_options.Inter);
    subplot(4,1,4);
    snr(x_after_adc(1:length(x_after_adc)), sim_options.Fs/sim_options.Inter);

    snr_in_id(num) = snr(x_to_subadc, sim_options.Fs/sim_options.Inter);
    snr_input(num) = snr(x_after_subadc, sim_options.Fs/sim_options.Inter);
    snr_output(num) = snr(x_after_adc, sim_options.Fs/sim_options.Inter);

    sfdr_in_id(num) = sfdr(x_to_subadc, sim_options.Fs/sim_options.Inter);
    sfdr_input(num) = sfdr(x_after_subadc, sim_options.Fs/sim_options.Inter);
    sfdr_output(num) = sfdr(x_after_adc, sim_options.Fs/sim_options.Inter);
    norm_freq(num) = sim_options.freq/(sim_options.Fs/sim_options.Inter/sim_options.M);
end

    %% Measurements2
    figure(7);
    subplot(2,1,1)
    plot(norm_freq, snr_in_id, '-o', norm_freq, snr_input, '-o', norm_freq, snr_output, '-o');
    title('SNR')
    xlabel('Нормированная частота') 
    ylabel('SNR (dB)') 
    legend('до калибровки без искажений', 'до калибровки с искажениями', 'после калибровки')
    subplot(2,1,2)
    plot(norm_freq, sfdr_in_id, '-o', norm_freq, sfdr_input, '-o', norm_freq, sfdr_output, '-o');
    title('SFDR (dB)')
    xlabel({'Нормированная частота fнорм = f/(Fs/M)','Fs - частота дискретизации всего TI-ADC, М - количество каналов'}) 
    ylabel('SFDR (dB)') 
    legend('до калибровки без искажений', 'до калибровки с искажениями','после калибровки')

    x4 = xline(0.42, '--', 'Интервал из статьи 1-ой зоны Найквиста')
    x4.LabelHorizontalAlignment = 'center'
    x4.LabelVerticalAlignment = 'middle';
    x2 = xline(0.55, '--', 'Интервал из статьи начало 2-ой зоны Найквиста')
    x2.LabelHorizontalAlignment = 'center'
    x2.LabelVerticalAlignment = 'middle';
    x3 = xline(0.92, '--', 'Интервал из статьи конец 2-ой зоны Найквиста')
    x3.LabelHorizontalAlignment = 'center'
    x3.LabelVerticalAlignment = 'middle';
    y2 = yline(79,'--', 'Нижняя граница SFDR (dB)')
    y2.LabelHorizontalAlignment = 'left'


stop_time = clock;
elapsed_time = etime(stop_time,start_time);

fprintf('Simulation duration: %g seconds\n',elapsed_time);