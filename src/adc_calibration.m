
% 1) Hu.M, Yi.P, (2022), Digital Calibration for Gain, Time Skew, and Bandwidth Mismatch 
%    in Under-Sampling Time-Interleaved System
% 2) Джиган В.И, Адаптивные фильтры
% 3) Айфичер Э, Джервис Б, Цифровая обработка сигналов. Практический подход
% 4) Behrouz Farhang-Boroujeny, Adaptive Filters Theory and Applications 

function [sig_adc, x_after_adc, x_after_adc_int] = adc_calibration(sim_options, adc_input, adc_input_int, s_to_subadc, s_after_subadc, Z)
    %% Calibration algorithm 1 (Fractional delays)

	% Fractional delays of ADC0 signal
	[yri_cut, yri_cut_int, yri_cut1, sig_adc] = fractional_delays(adc_input, adc_input_int, sim_options.M, sim_options.N, Z);

    %% Calibration algorithm 2 (Least Mean Squares)

    % adc_input = double(adc_input_int)*2^-11;
    % yri_cut = double(yri_cut_int)*2^-11; 

    adc_input = double(adc_input_int);
    yri_cut = double(yri_cut_int);

    [y_array, y_array_int] = least_mean_squares(adc_input, adc_input_int, yri_cut, yri_cut_int, sim_options.M, sim_options.N1, sim_options.Width);

    % create main signal after LS algorithm (switch after sub-adc)
    x_after_adc = zeros(length(y_array)*sim_options.M,1);
    x_after_adc_int = zeros(length(y_array_int)*sim_options.M,1);

    for i = 1:sim_options.M
        if i == 1
            % x_after_adc(i:sim_options.M:end) = yri_cut1(1:length(y_array),1);
            x_after_adc(i:sim_options.M:end) = double(yri_cut_int(1:length(y_array),1));
            x_after_adc_int(i:sim_options.M:end) = double(yri_cut_int(1:length(y_array),1));
        else
            x_after_adc(i:sim_options.M:end) = y_array(:,i-1);
            x_after_adc_int(i:sim_options.M:end) = double(y_array_int(:,i-1)) * 2^-sim_options.Width;
        end
    end

    figure(3);
    subplot(2,1,1)
    plot([x_after_adc_int]);
    title('Выход адаптивного фильтра int')
    xlabel('Номер отсчета') 
    ylabel('Амплитуда') 

    subplot(2,1,2)
    plot([x_after_adc]);
    title('Выход адаптивного фильтра double')
    xlabel('Номер отсчета') 
    ylabel('Амплитуда');
    % legend('double', 'double', 'Исходный сигнал с ошибками')

    % figure(4);
    % subplot(4,1,1);
    % snr(s_to_subadc(1:length(s_to_subadc)), sim_options.Fs/sim_options.Inter);
    % subplot(4,1,2);
    % snr(s_after_subadc(1:length(s_after_subadc)), sim_options.Fs/sim_options.Inter);
    % subplot(4,1,3);
    % snr(x_after_adc(1:length(x_after_adc)), sim_options.Fs/sim_options.Inter);
    % subplot(4,1,4);
    % snr(x_after_adc_int(1:length(x_after_adc_int)), sim_options.Fs/sim_options.Inter);

end
