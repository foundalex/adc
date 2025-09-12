
% 1) Hu.M, Yi.P, (2022), Digital Calibration for Gain, Time Skew, and Bandwidth Mismatch 
%    in Under-Sampling Time-Interleaved System
% 2) Джиган В.И, Адаптивные фильтры
% 3) Айфичер Э, Джервис Б, Цифровая обработка сигналов. Практический подход
% 4) Behrouz Farhang-Boroujeny, Adaptive Filters Theory and Applications 

function [sig_adc, x_after_adc] = adc_calibration(sim_options, adc_input, adc_input_int, s_after_subadc)
    %% Calibration algorithm 1.1 (Fractional delays)

	% Fractional delays of ADC0 signal
	[yri_cut, yri_cut_int, yri_cut1, sig_adc] = fractional_delays(adc_input, adc_input_int, sim_options.M, sim_options.N, sim_options.Z);

    %% Calibration algorithm 1.2 (Least Mean Squares)
    [y_array, y_array_int, error_out] = least_mean_squares(adc_input, adc_input_int, yri_cut, yri_cut_int, sim_options.M, 5);

    % create main signal after LS algorithm (switch after sub-adc)
    x_after_adc = zeros(length(y_array)*sim_options.M,1);
    x_after_adc_int = zeros(length(y_array_int)*sim_options.M,1);

    for i = 1:sim_options.M
        if i == 1
            x_after_adc(i:sim_options.M:end) = yri_cut1(1:length(y_array),1);
            x_after_adc_int(i:sim_options.M:end) = double(yri_cut_int(1:length(y_array),1))*2^-10;
        else
            x_after_adc(i:sim_options.M:end) = y_array(:,i-1);
            x_after_adc_int(i:sim_options.M:end) = double(y_array_int(:,i-1))*2^-11;
        end
    end

        % for i = 1:sim_options.M
        %     x_after_adc(i:sim_options.M:end) = yri_cut1(:,i);
        % end
    % end

    figure(3);
    plot([x_after_adc(1:250), x_after_adc_int(1:250)]);

    figure(4);
    subplot(3,1,1);
    snr(x_after_adc(1:length(x_after_adc)), sim_options.Fs/sim_options.Inter);
    subplot(3,1,2);
    snr(s_after_subadc(1:length(s_after_subadc)), sim_options.Fs/sim_options.Inter);
    subplot(3,1,3);
    snr(x_after_adc_int(1:length(x_after_adc_int)), sim_options.Fs/sim_options.Inter);
end
