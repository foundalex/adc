
% 1) Hu.M, Yi.P, (2022), Digital Calibration for Gain, Time Skew, and Bandwidth Mismatch 
%    in Under-Sampling Time-Interleaved System
% 2) Джиган В.И, Адаптивные фильтры
% 3) Айфичер Э, Джервис Б, Цифровая обработка сигналов. Практический подход
% 4) Behrouz Farhang-Boroujeny, Adaptive Filters Theory and Applications 

function [x_after_adc, x_after_adc_int, snr_s, fractional_mult, fractional_sum, fractional_width_total_mult, fractional_width_total_sum, ...  
    hilbert_mult, hilbert_sum, hilbert_width_total_mult, hilbert_width_total_sum, ...
    ... % Determinant
    Det2x2_mult1_abs, Det2x2_mult2_abs, Det2x2_sum_abs, ...
    Mult_DetM_2x2_1_abs_max, ...
	Mult_DetM_2x2_2_abs_max, ...
	Mult_DetM_2x2_3_abs_max, ...
	Mult_DetM_2x2_4_abs_max, ...
	Mult_DetM_2x2_5_abs_max, ...
	Mult_DetM_2x2_6_abs_max, ...
	Mult_DetM_2x2_7_abs_max, ...
	Mult_DetM_2x2_8_abs_max, ...
	Mult_DetM_2x2_9_abs_max, ...
	Mult_DetM_2x2_10_abs_max, ...
	... %% сумматоры определителя 3х3
	DetM_3x3_n_11_int_sum1_abs_max, ...
	DetM_3x3_n_11_int_abs_max, ...
	DetM_3x3_n_12_int_sum1_abs_max, ...
	DetM_3x3_n_12_int_abs_max, ...
	DetM_3x3_n_13_int_sum1_abs_max, ...
	DetM_3x3_n_13_int_abs_max, ...
	DetM_3x3_n_14_int_sum1_abs_max, ...
	DetM_3x3_n_14_int_abs_max, ...
	DetM_3x3_n_22_int_sum1_abs_max, ...
	DetM_3x3_n_22_int_abs_max, ...
	DetM_3x3_n_23_int_sum1_abs_max, ...
	DetM_3x3_n_23_int_abs_max, ...
	DetM_3x3_n_24_int_sum1_abs_max, ...
	DetM_3x3_n_24_int_abs_max, ...
	DetM_3x3_n_33_int_sum1_abs_max, ...
	DetM_3x3_n_33_int_abs_max, ...
	DetM_3x3_n_34_int_sum1_abs_max, ...
	DetM_3x3_n_34_int_abs_max, ...
	DetM_3x3_n_44_int_sum1_abs_max, ...
	DetM_3x3_n_44_int_abs_max ...
] = adc_calibration(sim_options, adc_input, s_to_subadc, s_after_subadc)

    %% Calibration algorithm 1 (Fractional delays)
    n = (0:1:sim_options.N-1);
    Nbp = floor(sim_options.Z/2); % стр 7. (24)
    del_proc = ((sim_options.N-1)/2);
    nn = 1:length(adc_input(:,1));
    w_blackman = 0.42 - 0.5 * cos(2*pi*n/(sim_options.N-1)) + 0.08 * cos(4*pi*n/(sim_options.N-1)); % Blackman window

    %% Fractional filter coeff
    %%
    delay_adc = (1/sim_options.M:1/sim_options.M:1); % (стр.6,(16)), создаем массив на различные значения задержек
    D = del_proc - delay_adc; % delay (N-1)/2 - d = causal filter
    hri_m = sinc(n'- D); % shift impulse response on D = Dint - d for fractional delay filter
    w_blackman_fractional = 0.42 - 0.5 * cos(2*pi*(n'+ delay_adc)/(sim_options.N-1)) + 0.08 * cos(4*pi*(n'+ delay_adc)/(sim_options.N-1)); % shift Blackman window
    hri_w = hri_m .* w_blackman_fractional; 

    
    coeff_frac_int = cast((hri_w*2^(sim_options.fractional_coeff_width-1)), sim_options.int_size);
    
    % figure(3)
    % subplot(2,1,1)
    % plot(double(coeff_frac_int))
    % title('Импульсная характеристика фильтра Гилберта')
    % xlabel('Номер отсчета') 
    % ylabel('Амплитуда') 
    % subplot(2,1,2)
    % plot(width)
    % title('Разрядность коэффициентов')
    % xlabel('Номер коэффициента') 
    % ylabel('Необходимое количество бит') 

    % [y, f] = freqz(hri_w(:,1), 1,1024, 'whole', 1000000000);
    % for k = 1:sim_options.M-1
    %     hrim_fi_test = fi(hri_w(:,k), 1,fractional_width(1),fractional_width(1)-1);
    %     hrim_int_test = int32(hrim_fi_test * 2^(fractional_width(1)-1));
    %     [y1(:,k), f1(:,k)] = freqz(double(hrim_int_test)*2^-(fractional_width(1)-1),1,1024, 'whole', 1000000000);
    % end
    % 
    % figure(2);
    % plot(f, abs(y), f, abs(y1(:,1)), f, abs(y1(:,2)), f, abs(y1(:,3)));
    % title('Влияние разрядностей коэффициентов на АЧХ фильтра дробной задержки')
    % xlabel('Частота') 
    % ylabel('Коэффициент передачи') 
    % legend({'double','16 бит', '19 бит', '21 бит'},'Location','northeast')

    %% Hilbert filter coeff
    %%
    hh = (2./((n-del_proc)*pi)).*(sin(((n-del_proc)*pi)./2)).^2;
    hh(1) = 0;
    hh(37) = 0;
    hh_m = (hh .* w_blackman).';

    % Negative Symmetric coefficients
    hilbert_coeff_int = cast(hh_m*2^(sim_options.hilbert_coeff_width-1), sim_options.int_size);

    % [y, f] = freqz(double(hilbert_coeff_int)*2^-(hilbert_width-1), 1,1024, 'whole', 1000000000);
    % figure(2);
    % plot(f, abs(y));
    % 
    % figure(3)
    % subplot(2,1,1)
    % plot(double(hh_m_int))
    % title('Импульсная характеристика фильтра Гилберта')
    % xlabel('Номер отсчета') 
    % ylabel('Амплитуда') 
    % subplot(2,1,2)
    % plot(width)
    % title('Разрядность коэффициентов')
    % xlabel('Номер коэффициента') 
    % ylabel('Необходимое количество бит') 

    %% zones Nyquist
    nn1 = nn' + delay_adc;
    a1 = 2*pi*nn1*Nbp;
    cosi = cos(a1);
    sini = sin(a1);
    % sini = sin(a1(1:end-del_proc,:));

    %%
    yri_cut = zeros(length(adc_input(:,1)),sim_options.M);
    yri_cut_int = cast(zeros(length(adc_input(:,1)),sim_options.M), sim_options.int_size);

    fractional_mult = cast(zeros(sim_options.N, sim_options.M-1), sim_options.int_size); 
    fractional_sum = cast(zeros(sim_options.N-1, sim_options.M-1), sim_options.int_size); 
    fractional_width_total_mult = cast(zeros(sim_options.N, sim_options.M-1), sim_options.int_size); 
    fractional_width_total_sum = cast(zeros(sim_options.N-1, sim_options.M-1), sim_options.int_size); 

    hilbert_mult = cast(zeros(sim_options.N, sim_options.M-1), sim_options.int_size); 
    hilbert_sum = cast(zeros(sim_options.N-1, sim_options.M-1), sim_options.int_size); 
    hilbert_width_total_mult = cast(zeros(sim_options.N, sim_options.M-1), sim_options.int_size); 
    hilbert_width_total_sum = cast(zeros(sim_options.N-1, sim_options.M-1), sim_options.int_size); 

	%% Fractional delays of ADC0 signal
    % for i = 1:sim_options.M-1

    %     fractional_mult_file = ['src/width_txt/Width_multiplier_Fractional_filter_' num2str(i) '.txt'];
    %     fractional_sum_file = ['src/width_txt/Width_adder_Fractional_filter_' num2str(i) '.txt'];
    % 
    %     hilbert_mult_file = ['src/width_txt/Width_multiplier_Hilbert_filter_' num2str(i) '.txt'];
    %     hilbert_sum_file = ['src/width_txt/Width_adder_Hilbert_filter_' num2str(i) '.txt'];
    % 
    %     yri = filter(hri_w(:,i), 1, adc_input(:,1)); % filter (стр.6 (15))
    %     ymi = filter(hh_m.', 1, yri);
    % 
    %     % y(n) = x(n)*(k1*2^N)+x(n-1)*(k2*2^N)+x(n-3)*(k3*2^N)
    %     % N = 18
    %     [y_fractional_outInt, fractional_mult(:,i), fractional_sum(:,i), fractional_width_total_mult(:,i), fractional_width_total_sum(:,i) ...
    %         ] = fir_filter(coeff_frac_int(:,i), adc_input(:,1), sim_options.N, ...
    %     fractional_mult_file, fractional_sum_file, sim_options.width_fractional, sim_options); % (стр.6 (15)) 
    % 
    %     %% Hilbert
    %     %%
    %     % y(n) = x(n)*(k1*2^N)+x(n-1)*(k2*2^N)+x(n-3)*(k3*2^N)
    %     % N = 15
    %     [ymi_HilbertInt, hilbert_mult(:,i), hilbert_sum(:,i), hilbert_width_total_mult(:,i), hilbert_width_total_sum(:,i) ...
    %         ] = fir_filter(hilbert_coeff_int, y_fractional_outInt, sim_options.N, ...
    %      hilbert_mult_file, hilbert_sum_file, sim_options.width_hilbert, sim_options); % (стр.6 (15)) );
    % 
    %     %%
    % 
    %     snr_fractional_out_double = snr(yri, 1000000000);
    %     snr_fractional_out_int = snr(double(y_fractional_outInt)*2^-(sim_options.fractional_coeff_width-1), 1000000000);
    % 
    %     if (snr_fractional_out_double - snr_fractional_out_int) > 0.1
    %         disp('SNR fractional out different!')
    %         disp([sim_options.SNR, sim_options.freq])
    %     end
    % 
    %     y_fractional_outInt_double = double(y_fractional_outInt)*2^-(sim_options.fractional_coeff_width-1);
    %     relative_error_fractional = yri./y_fractional_outInt_double;
    % 
    %     sim_options.divide_remainder = (sim_options.fractional_coeff_width-1) + (sim_options.hilbert_coeff_width-1);
    % 
    %     ymi_HilbertInt_double = double(ymi_HilbertInt)*2^-(sim_options.divide_remainder);
    %     relative_error_hilbert = ymi./ymi_HilbertInt_double;
    % 
    %     snr_hilbert_out_double = snr(ymi, 1000000000);
    %     snr_hilbert_out_int = snr(double(ymi_HilbertInt)*2^-30, 1000000000);
    % 
    % 
    %     if (snr_hilbert_out_double - snr_hilbert_out_int) > 0.1
    %         disp('SNR hilbert out different!')
    %         disp([sim_options.SNR, sim_options.freq])
    %     end
    % 
    %     figure(4);
    %     subplot(4,1,1)
    %     plot([yri(1:500), y_fractional_outInt_double(1:500)]);
    %     subplot(4,1,2);
    %     snr(yri, 1000000000);
    %     subplot(4,1,3);
    %     snr(y_fractional_outInt_double, 1000000000);
    % 
    %     subplot(4,1,4);
    %     plot(relative_error_fractional);
    %     title('Относительная ошибка выходного сигнала фильтра дробной задержки между double и integer')
    %     xlabel('Номер отсчета') 
    %     ylabel('Значение ошибки') 
    %     x4 = xline(37, '--', 'Переходной процесс фильтра')
    %     x4.LabelHorizontalAlignment = 'center'
    %     x4.LabelVerticalAlignment = 'middle';
    % 
    %     figure(5);
    %     subplot(4,1,1)
    %     plot([ymi(1:500), double(ymi_HilbertInt(1:500))*2^-sim_options.divide_remainder]);
    %     subplot(4,1,2);
    %     snr(ymi, 1000000000);
    %     subplot(4,1,3);
    %     snr(double(ymi_HilbertInt)*2^-sim_options.divide_remainder, 1000000000);
    %     subplot(4,1,4);
    %     plot(relative_error_hilbert);
    %     title('Относительная ошибка выходного сигнала фильтра Гилберта между double и integer')
    %     xlabel('Номер отсчета') 
    %     ylabel('Значение ошибки') 
    %     x4 = xline(73, '--', 'Переходной процесс фильтра')
    %     x4.LabelHorizontalAlignment = 'center'
    %     x4.LabelVerticalAlignment = 'middle';
    % 
    %     %%
    %     [yric(:,i), yric_int(:,i)] = single_sideband(yri, y_fractional_outInt, ymi, ymi_HilbertInt, cosi(:,i), sini(:,i), i, del_proc, sim_options);
    % 
    % end
    % 
    % 
    % %% test signal
    % yri_cut = zeros(length(adc_input(1:end-del_proc,1)),sim_options.M);
    % yri_cut_int = cast(zeros(length(adc_input(1:end-del_proc,1)),sim_options.M), sim_options.int_size);
    % 
    % 
    % for i = 1:sim_options.M
    %     if i == 1
    %          yri_cut(:,1) = adc_input(1:end-del_proc,1);
    %          yri_cut_int(:,1) = adc_input(1:end-del_proc,1);
    %     else
    %         yri_cut(:,i) = yric(del_proc+1:end,i-1);
    %         yri_cut_int(:,i) = yric_int(del_proc+1:end,i-1);
    %     end
    % end
    % 
    % yri_cut(end-del_proc:end,:) = [];
    % yri_cut_int(end-del_proc:end,:) = [];
    % 
    % 
    % sig_adc = zeros(sim_options.M*length(yri_cut(:,1)),1);
    % sig_adc_int = zeros(sim_options.M*length(yri_cut_int(:,1)),1);
    % 
    % yri_cut1(:,1) = double(yri_cut_int(:,1));
    % yri_cut1(:,2) = double(yri_cut_int(:,2))*2^-(sim_options.fractional_coeff_width-1+sim_options.hilbert_coeff_width-1);
    % yri_cut1(:,3) = double(yri_cut_int(:,3))*2^-(sim_options.fractional_coeff_width-1);
    % yri_cut1(:,4) = double(yri_cut_int(:,4))*2^-(sim_options.fractional_coeff_width-1+sim_options.hilbert_coeff_width-1);
    % % 
    % % 
	% for i = 1:sim_options.M
    %     sig_adc(i:sim_options.M:end) = yri_cut(:,i);
    %     sig_adc_int(i:sim_options.M:end) = yri_cut1(:,i);
    % end




    
    % 
    % save (sprintf(num2str(clock) + ".mat"));
    % load ('2025             11              5             11              8         52.084.mat'); % 888 MHz 70 SNR

    % figure(6);
    % plot([s_to_subadc(1:500), sig_adc(1:500), sig_adc_int(1:500)]);
    % figure(10);
    % subplot(3,1,1)
    % snr(s_to_subadc, sim_options.Fs/sim_options.Inter);
    % subplot(3,1,2);
    % snr(sig_adc, sim_options.Fs/sim_options.Inter);
    % subplot(3,1,3);
    % snr(sig_adc_int, sim_options.Fs/sim_options.Inter);
    %%

    snr_s = 0;

    %% Calibration algorithm 2 (Least Mean Squares)

    % for i = 1:sim_options.M-1
    %     [y_array(:,i), y_array_int(:,i), DetM_2x2_array(:,i), DetM_2x2_array_int(:,i), Det2x2_mult1_abs(:,i), Det2x2_mult2_abs(:,i), Det2x2_sum_abs(:,i), ...
    %         Mult_DetM_2x2_1_abs_max(:,i), ...
	%         Mult_DetM_2x2_2_abs_max(:,i), ...
	%         Mult_DetM_2x2_3_abs_max(:,i), ...
	%         Mult_DetM_2x2_4_abs_max(:,i), ...
	%         Mult_DetM_2x2_5_abs_max(:,i), ...
	%         Mult_DetM_2x2_6_abs_max(:,i), ...
	%         Mult_DetM_2x2_7_abs_max(:,i), ...
	%         Mult_DetM_2x2_8_abs_max(:,i), ...
	%         Mult_DetM_2x2_9_abs_max(:,i), ...
	%         Mult_DetM_2x2_10_abs_max(:,i), ...
	% 		... %% сумматоры определителя 3х3
	% 		DetM_3x3_n_11_int_sum1_abs_max(:,i), ...
	% 		DetM_3x3_n_11_int_abs_max(:,i), ...
	% 		DetM_3x3_n_12_int_sum1_abs_max(:,i), ...
	% 		DetM_3x3_n_12_int_abs_max(:,i), ...
	% 		DetM_3x3_n_13_int_sum1_abs_max(:,i), ...
	% 		DetM_3x3_n_13_int_abs_max(:,i), ...
	% 		DetM_3x3_n_14_int_sum1_abs_max(:,i), ...
	% 		DetM_3x3_n_14_int_abs_max(:,i), ...
	% 		DetM_3x3_n_22_int_sum1_abs_max(:,i), ...
	% 		DetM_3x3_n_22_int_abs_max(:,i), ...
	% 		DetM_3x3_n_23_int_sum1_abs_max(:,i), ...
	% 		DetM_3x3_n_23_int_abs_max(:,i), ...
	% 		DetM_3x3_n_24_int_sum1_abs_max(:,i), ...
	% 		DetM_3x3_n_24_int_abs_max(:,i), ...
	% 		DetM_3x3_n_33_int_sum1_abs_max(:,i), ...
	% 		DetM_3x3_n_33_int_abs_max(:,i), ...
	% 		DetM_3x3_n_34_int_sum1_abs_max(:,i), ...
	% 		DetM_3x3_n_34_int_abs_max(:,i), ...
	% 		DetM_3x3_n_44_int_sum1_abs_max(:,i), ...
	% 		DetM_3x3_n_44_int_abs_max(:,i) ...
    %     ] = least_mean_squares(double(adc_input(:,i+1)), adc_input(:,i+1), yri_cut(:,i+1), yri_cut_int(:,i+1), sim_options);
    % 
    % 
    %     %%  
    % 
    %     figure(13);
    % 
    %     subplot(2,1,1);
    %     plot(double(yri_cut_int(:,1)));
    %     title('Выход референсного канала АЦП')
    %     xlabel('Номер отсчета') 
    %     ylabel('Амплитуда сигнала') 
    % 
    %     subplot(2,1,2);
    %     plot([y_array(:,i)]);
    %     title('Выход адаптивного фильтра double для одного канала АЦП')
    %     xlabel('Номер отсчета') 
    %     ylabel('Амплитуда сигнала') 
    % 
    % end
    



    % save (sprintf(num2str(clock) + ".mat"));
    load ('2025             11              5             11             34         46.136.mat'); % 888 MHz 70 SNR









    % % create main signal after LS algorithm (switch after sub-adc)
    % x_after_adc = zeros(length(y_array)*sim_options.M,1);
    % x_after_adc_int = zeros(length(y_array_int)*sim_options.M,1);
    % 
    % for i = 1:sim_options.M
    %     if i == 1
    %         x_after_adc(i:sim_options.M:end) = yri_cut1(1:length(y_array),1);
    %         % x_after_adc(i:sim_options.M:end) = double(yri_cut_int(1:length(y_array),1));
    %         x_after_adc_int(i:sim_options.M:end) = double(yri_cut_int(1:length(y_array),1));
    %     else
    %         x_after_adc(i:sim_options.M:end) = y_array(:,i-1);
    %         x_after_adc_int(i:sim_options.M:end) = double(y_array_int(:,i-1)) * 2^-32;
    %     end
    % end
    % % 
    % figure(3);
    % subplot(2,1,1)
    % plot([x_after_adc_int]);
    % title('Выход адаптивного фильтра int')
    % xlabel('Номер отсчета') 
    % ylabel('Амплитуда') 
    % 
    % subplot(2,1,2)
    % plot([x_after_adc]);
    % title('Выход адаптивного фильтра double')
    % xlabel('Номер отсчета') 
    % ylabel('Амплитуда');
    % 
    % 
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
