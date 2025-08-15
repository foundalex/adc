
% 1) Hu.M, Yi.P, (2022), Digital Calibration for Gain, Time Skew, and Bandwidth Mismatch 
%    in Under-Sampling Time-Interleaved System
% 2) Джиган В.И, Адаптивные фильтры
% 3) Айфичер Э, Джервис Б, Цифровая обработка сигналов. Практический подход
% 4) Behrouz Farhang-Boroujeny, Adaptive Filters Theory and Applications 

function [sig_adc, x_after_adc, error_out] = adc_calibration(sim_options, adc_input)
    %% Calibration algorithm 1.1 (Fractional delays)

	% Fractional delays of ADC0 signal
	[yri_cut] = fractional_delays(adc_input, sim_options.M, sim_options.N, sim_options.Z);

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

	% figure(2);
	% subplot(3,1,1);
	% snr(x_to_subadc(1:length(sig_adc)), Fs/(sim_options.Inter));
	% subplot(3,1,2);
	% snr(x_after_subadc(1:length(sig_adc)), Fs/(sim_options.Inter));
	% subplot(3,1,3);
	% snr(sig_adc(1000:length(sig_adc)), Fs/(sim_options.Inter));

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
    [y_array, error_out] = least_mean_squares(adc_input, yri_cut, sim_options.M, sim_options.N);

    % create main signal after LS algorithm (switch after sub-adc)
    x_after_adc = zeros(length(y_array)*sim_options.M,1);
    for i = 1:sim_options.M
        if i == 1
            x_after_adc(i:sim_options.M:end) = yri_cut(1:length(y_array));
        else
            x_after_adc(i:sim_options.M:end) = y_array(:,i-1);
        end
    end

end
