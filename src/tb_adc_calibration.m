% 1) Джиган В.И, Адаптивные фильтры
% 2) Айфичер Э, Джервис Б, Цифровая обработка сигналов. Практический подход
% 3) Hu.M, Yi.P, (2022), Digital Calibration for Gain, Time Skew, and Bandwidth Mismatch 
%    in Under-Sampling Time-Interleaved System
% 4) Behrouz Farhang-Boroujeny, Adaptive Filters Theory and Applications 


% close all;
% clear all;
% clc;

%% Parameters;
global sim_consts;

% M = 8;                        % Num of sub ADC
N = 73;                         % Number taps filters 
StopTime = 0.00001;             % seconds
Inter = 10;                     % oversampling factor
Fs = 1000000000 * Inter * sim_options.M;    % Fs all ADC system
dt = 1/Fs;                      % seconds per sample
t = (0:dt:StopTime-dt)';        % seconds
% sim_options.freq = 105000000;               % initial tone frequency
num_cycles = 1;                 % number of cycles testing algorithm
%% Error parameters
% offset_error_array = [0.5 0.4 0.3 0.5 0.4 0.3 0.2];        % empty offset
time_skew_array = [0.5 0.2 0.1 0.4 0.3 0.1 0.2]; % time skew error array of sub-adc (ADC2-ADC8) (fc/time_skew_array(i))
gain_error_array = [1.4 1.3 1.2 1.1 1.2 1.3 1.1]; % gain error array of sub-adc (ADC2-ADC8)
MODEL_ERROR = true;
%% tb
% freq_array = [225000000 405000000];
% freq = gpuArray(freq);
% isgpuarray(freq)

for num = 1:num_cycles

    step = 10000000;                 % step frequency input signal
    % freq = sim_options.freq + step;              % frequency of fundamental tone
    Z = ceil(sim_options.freq/(Fs/Inter/2/sim_options.M));   % Nyquist zone

    % Create main signal with noise
    s = 0.75*cos(2*pi*sim_options.freq*t);
    noise = awgn(s,sim_options.SNR);
    s = s + noise;

    % oversampled signal transfer to sub-adc
    adc_input(:,1) = s(1:sim_options.M*Inter:end);
    for i = 1:sim_options.M-1
        adc_input(:,i+1) = s(i*Inter+1:sim_options.M*Inter:end);
        indexx(i) = (i*Inter+1);
    end


	% исходный сигнал до искажений
	x_to_subadc = zeros(sim_options.M*length(adc_input(:,1)),1);
	for i = 1:sim_options.M
		x_to_subadc(i:sim_options.M:end) = adc_input(:,i); 
	end

    %% Add time skew error, gain error
    if MODEL_ERROR == true
        num_adc = sim_options.M;
        % time skew model
        for i = 1:sim_options.M-1
            adc_input_skew = time_skew_func(time_skew_array(i), s, indexx(i), Inter, num_adc); 
            adc_input(:,i+1) = adc_input_skew;
        end

	    % % model offset error
        % for i = 1:M-1
        %     adc_input(:,i+1) = adc_input(:,i+1) + offset_error_array(i);
        % end


        % model gain error
        for i = 1:sim_options.M-1
            adc_input(:,i+1) = adc_input(:,i+1) * gain_error_array(i);
        end
    end

    % Main signal with gain error and time skew
	x_after_subadc = zeros(sim_options.M*length(adc_input(:,1)),1);
	for i = 1:sim_options.M
	    x_after_subadc(i:sim_options.M:end) = adc_input(:,i); 
    end
    %% Calibration algorithm 1.1 (Fractional delays)
	% Hu.M, Yi.P, (2022), Digital Calibration for Gain, Time Skew, and Bandwidth Mismatch 
    % in Under-Sampling Time-Interleaved System

	% Fractional delays of ADC0 signal
	[yri_cut] = fractional_delays(adc_input, sim_options.M, N, Z);

    %% test signal after fractional delay filters
	sig_adc = zeros(sim_options.M*length(yri_cut(:,1)),1);
	for i = 1:sim_options.M
		sig_adc(i:sim_options.M:end) = yri_cut(:,i);
    end

	% figure(11);
	% plot([x_after_subadc(1:500), sig_adc(1:500)]);

	% figure(1);
	% subplot(3,1,1);
	% sfdr(x_to_subadc(1:length(sig_adc)), Fs/(Inter));
	% subplot(3,1,2);
	% sfdr(x_after_subadc(1:length(sig_adc)), Fs/(Inter));
	% subplot(3,1,3);
	% sfdr(sig_adc(1000:length(sig_adc)), Fs/(Inter));
	% % 
	figure(2);
	subplot(3,1,1);
	snr(x_to_subadc(1:length(sig_adc)), Fs/(Inter));
	subplot(3,1,2);
	snr(x_after_subadc(1:length(sig_adc)), Fs/(Inter));
	subplot(3,1,3);
	snr(sig_adc(1000:length(sig_adc)), Fs/(Inter));

	% spectrumScope = spectrumAnalyzer(SampleRate=Fs, ...            
	%             AveragingMethod='exponential',ForgettingFactor=0, ...
	%             YLimits=[-30 10],ShowLegend=true, Method='Welch');
    % 
	% spectrumScope.WindowLength = 2048;
	% spectrumScope.FrequencyResolutionMethod = "window-length";
	% spectrumScope.PlotAsTwoSidedSpectrum=false;
	% spectrumScope.DistortionMeasurements.Enabled = true;
    % 
	% spectrumScope([sig_adc(4096:19000)]);

    %% Calibration algorithm 1.2 (Least Mean Squares)
    [y_array, error_out] = least_mean_squares(adc_input, yri_cut, sim_options.M, N);

    % create main signal after LS algorithm (switch after sub-adc)
    x_after_adc = zeros(length(y_array)*sim_options.M,1);
    for i = 1:sim_options.M
        if i == 1
            x_after_adc(i:sim_options.M:end) = yri_cut(1:length(y_array));
        else
            x_after_adc(i:sim_options.M:end) = y_array(:,i-1);
        end
    end

    %% Measurements
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
    sfdr(x_to_subadc(1:length(x_after_adc)), Fs/Inter);
    subplot(4,1,2);
    sfdr(x_after_subadc(1:length(x_after_adc)), Fs/Inter);
    subplot(4,1,3);
    sfdr(sig_adc(1:length(x_after_adc)), Fs/Inter);
    subplot(4,1,4);
    sfdr(x_after_adc(1:length(x_after_adc)), Fs/Inter);

    figure(6);
    subplot(4,1,1);
    snr(x_to_subadc(1:length(x_after_adc)), Fs/Inter);
    subplot(4,1,2);
    snr(x_after_subadc(1:length(x_after_adc)), Fs/Inter);
    subplot(4,1,3);
    snr(sig_adc(1:length(x_after_adc)), Fs/Inter);
    subplot(4,1,4);
    snr(x_after_adc(1:length(x_after_adc)), Fs/Inter);

    snr_in_id(num) = snr(x_to_subadc, Fs/Inter);
    snr_input(num) = snr(x_after_subadc, Fs/Inter);
    snr_output(num) = snr(x_after_adc, Fs/Inter);

    sfdr_in_id(num) = sfdr(x_to_subadc, Fs/Inter);
    sfdr_input(num) = sfdr(x_after_subadc, Fs/Inter);
    sfdr_output(num) = sfdr(x_after_adc, Fs/Inter);
    norm_freq(num) = sim_options.freq/(Fs/(Inter)/sim_options.M);

end

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
%%
% function for model timing skew

function adc_input_skew = time_skew_func(time_skew, s, indexx, Inter, num_adc) 
    adc_input_skew = s(indexx + time_skew*Inter:num_adc*Inter:end);
end