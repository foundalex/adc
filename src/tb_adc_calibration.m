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
            
for num = 1:sim_options.num_cycles

    % Функция генерации сигналов для АЦП
    [s_to_subadc, adc_input, s_after_subadc] = gen_oversampled_signal(sim_options);
    % Функция калибровки АЦП
    [sig_adc, delta_tilda] = new_algorithm(adc_input, s_to_subadc, s_after_subadc, sim_options);

    %% 
    snr_in_good(num) = snr(double(s_to_subadc), sim_options.Fs/sim_options.Inter); % SNR выходного сигнала без ошибок
    snr_in_bad(num) = snr(double(s_after_subadc), sim_options.Fs/sim_options.Inter); % SNR выходного сигнала с ошибками
    snr_output(num) = snr(sig_adc, sim_options.Fs/sim_options.Inter); % SNR выходного сигнала после калибровки

    sfdr_in_good(num) = sfdr(double(s_to_subadc), sim_options.Fs/sim_options.Inter);  % SFDR выходного сигнала без ошибок
    sfdr_in_bad(num) = sfdr(double(s_after_subadc), sim_options.Fs/sim_options.Inter); % SFDR выходного сигнала с ошибками
    sfdr_output(num) = sfdr(sig_adc, sim_options.Fs/sim_options.Inter); % SFDR выходного сигнала после калибровки

    norm_freq(num) = sim_options.freq/(sim_options.Fs/sim_options.Inter); % нормированная частота
    num_array(:,num) = num;

    % freq
    sim_options.freq = sim_options.freq + sim_options.step; % frequency of fundamental tone
end

% Итоговый график SNR и SFDR каждой итерации алгоритма
figure(14);
subplot(2,1,1)
plot(norm_freq, snr_in_good, '-o', norm_freq, snr_in_bad, '-o', norm_freq, snr_output, '-o');
title('SNR');
xlabel('Нормированная частота'); 
ylabel('SNR (dB)'); 
legend({'Входной сигнал без ошибок', 'Входной сигнал с ошибками', 'Выходной сигнал'}, 'Location','northwest');
% 
subplot(2,1,2)
plot(norm_freq, sfdr_in_good, '-o', norm_freq, sfdr_in_bad, '-o', norm_freq, sfdr_output, '-o');
title('SFDR (dB)');
xlabel('Нормированная частота'); 
ylabel('SFDR (dB)'); 
legend({'Входной сигнал без ошибок', 'Входной сигнал с ошибками', 'Выходной сигнал'}, 'Location','northwest');

%%
stop_time = clock;
elapsed_time = etime(stop_time,start_time);

fprintf('Simulation duration: %g seconds\n',elapsed_time);

if (sim_options.enable_log == true)
    diary off;
end