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

%% Фильтр дробной задержки
width_mult = int8(zeros(sim_options.N,sim_options.M-1));
width_sum = int8(zeros(sim_options.N-1,sim_options.M-1));
fractional_mult_max = cast(zeros(sim_options.N,sim_options.M-1), sim_options.int_size);
fractional_sum_max = cast(zeros(sim_options.N,sim_options.M-1), sim_options.int_size);
fractional_total_width_mult_max = cast(zeros(sim_options.N,sim_options.M-1), sim_options.int_size);
fractional_total_width_sum_max = cast(zeros(sim_options.N,sim_options.M-1), sim_options.int_size);
y_fractional_outInt_abs_max_in_cycle = cast(zeros(sim_options.M-1,1), sim_options.int_size);
%% Фильтр Гилберта
width_mult_h = int8(zeros(sim_options.N, sim_options.M-1));
width_sum_h = int8(zeros(sim_options.N-1, sim_options.M-1));
hilbert_mult_max = cast(zeros(sim_options.N,sim_options.M-1), sim_options.int_size);
hilbert_sum_max = cast(zeros(sim_options.N,sim_options.M-1), sim_options.int_size);
hilbert_total_width_mult_max = cast(zeros(sim_options.N,sim_options.M-1), sim_options.int_size);
hilbert_total_width_sum_max = cast(zeros(sim_options.N,sim_options.M-1), sim_options.int_size);
ymi_HilbertInt_abs_max_in_cycle = cast(zeros(sim_options.M-1,1), sim_options.int_size);
%%
num_filter_out = string((1:3)');
num_filter = string((1:sim_options.N)');
num_determinante = string((1:sim_options.num_det2x2*3)');
num_adaptive = string((1:sim_options.Size_matrix)');
num_divide = string((1:1)');

DetM_2x2_multiplier_total_abs_max_in_cycle 			= cast(zeros(sim_options.num_det2x2*3,sim_options.M-1), sim_options.type_2x2_det);
Det2x2_sum_abs_max_in_cycle 						= cast(zeros(sim_options.num_det2x2*3,sim_options.M-1), sim_options.type_2x2_det);
Mult_DetM_3x3_array_max_in_cycle 					= cast(zeros(sim_options.num_det2x2*3,sim_options.M-1), sim_options.type_3x3_det);
Mult_DetM_3x3_array_mult_total_width_max_in_cycle 	= cast(zeros(sim_options.num_det2x2*3,sim_options.M-1), sim_options.type_3x3_det);
DetM_3x3_int_pre_sum_array_max_in_cycle 			= cast(zeros(sim_options.num_det2x2*3,sim_options.M-1), sim_options.type_3x3_det);
DetM_3x3_int_pre_sum_width_total_max_in_cycle 		= cast(zeros(sim_options.num_det2x2*3,sim_options.M-1), sim_options.type_3x3_det);
DetM_3x3_int_sum_array_max_in_cycle 				= cast(zeros(sim_options.num_det2x2*3,sim_options.M-1), sim_options.type_3x3_det);
DetM_3x3_int_sum_width_total_max_in_cycle 			= cast(zeros(sim_options.num_det2x2*3,sim_options.M-1), sim_options.type_3x3_det);
% умножители определителя 4х4
DetM_4x4_int_mult_array_max_in_cycle 				= cast(zeros(sim_options.num_det2x2*3,sim_options.M-1), sim_options.type_4x4_det);
% разрядность умножителей определителя 4х4
DetM_4x4_int_mult_width_total_max_in_cycle 			= cast(zeros(sim_options.num_det2x2*3,sim_options.M-1), sim_options.type_4x4_det);
% пресумматоры определителя 4х4
DetM_4x4_int_pre_sum_array_max_in_cycle 			= cast(zeros(sim_options.num_det2x2*3,sim_options.M-1), sim_options.type_4x4_det);
% разрядность пресумматоров определителя 4х4
DetM_4x4_int_pre_sum_width_total_max_in_cycle 		= cast(zeros(sim_options.num_det2x2*3,sim_options.M-1), sim_options.type_4x4_det);
% сумматоры определителя 4х4
DetM_4x4_int_sum_array_max_in_cycle 				= cast(zeros(sim_options.num_det2x2*3,sim_options.M-1), sim_options.type_4x4_det);
% разрядность сумматоров определителя 4х4
DetM_4x4_int_sum_array_width_total_max_in_cycle		= cast(zeros(sim_options.num_det2x2*3,sim_options.M-1), sim_options.type_4x4_det);
% умножители определителя 5х5
DetM_5x5_int_mult_array_max_in_cycle 				= cast(zeros(sim_options.num_det2x2*3,sim_options.M-1), sim_options.type_4x4_det);
% разрядность умножителей 5x5
DetM_5x5_int_mult_array_width_total_max_in_cycle	= cast(zeros(sim_options.num_det2x2*3,sim_options.M-1), sim_options.type_5x5_det);
% пресумматоры1 определителя 5х5
DetM_5x5_int_pre_sum1_array_max_in_cycle 			= cast(zeros(sim_options.num_det2x2*3,sim_options.M-1), sim_options.type_5x5_det);
% разрядность пресумматоров1 определителя 5х5
DetM_5x5_int_pre_sum1_array_width_total_max_in_cycle = cast(zeros(sim_options.num_det2x2*3,sim_options.M-1), sim_options.type_5x5_det);
% пресумматор2 определителя 5х5
DetM_5x5_int_sum3_abs_max_in_cycle 					= cast(zeros(sim_options.num_det2x2*3,sim_options.M-1), sim_options.type_5x5_det);
% разрядность пресумматора2 определителя 5х5
DetM_5x5_int_sum3_width_total_max_in_cycle 			= cast(zeros(sim_options.num_det2x2*3,sim_options.M-1), sim_options.type_5x5_det);
% сумматор определителя 5х5
DetM_5x5_int_abs_max_in_cycle 						= cast(zeros(sim_options.num_det2x2*3,sim_options.M-1), sim_options.type_5x5_det);
% разрядность сумматора определителя 5х5
DetM_5x5_int_width_total_max_in_cycle 				= cast(zeros(sim_options.num_det2x2*3,sim_options.M-1), sim_options.type_5x5_det);
% начальный определитель
Det_x3_int_max_in_cycle                             = cast(zeros(1,sim_options.M-1), sim_options.type_5x5_det);
% выход делителя
Divide_max_in_cycle                                 = cast(zeros(1,sim_options.M-1), sim_options.type_divide_out);

% Адаптивный фильтр
Adaptive_filter_mult_array_max_in_cycle 	        = cast(zeros(sim_options.Size_matrix,sim_options.M-1), sim_options.type_mult_in_adaptive_filter);
Adaptive_filter_mult_total_width_in_cycle 	        = cast(zeros(sim_options.Size_matrix,sim_options.M-1), sim_options.type_mult_in_adaptive_filter);
Adaptive_filter_sum_array_max_in_cycle 	            = cast(zeros(sim_options.Size_matrix,sim_options.M-1), sim_options.type_add_in_adaptive_filter);
Adaptive_filter_sum_total_width_in_cycle 	        = cast(zeros(sim_options.Size_matrix,sim_options.M-1), sim_options.type_add_in_adaptive_filter);


for num = 1:sim_options.num_cycles

    % Функция генерации сигналов для АЦП
    [s_to_subadc, adc_input, s_after_subadc, sim_options.Z] = gen_oversampled_signal(sim_options);

    % Основная функция калибровки АЦП
    [x_after_adc, x_after_adc_double, x_after_adc_int, ...
        ... % Значения полосового фильтра ADC0
        golden_mult, golden_sum, golden_width_total_mult, golden_width_total_sum, y_golden_outInt_abs_max, ...  
        ... % Значения фильтров дробной задержки
        fractional_mult, fractional_sum, fractional_width_total_mult, fractional_width_total_sum, y_fractional_outInt_abs_max, ...  
        ... % Значения фильтров Гилберта
        hilbert_mult, hilbert_sum, hilbert_width_total_mult, hilbert_width_total_sum, ymi_HilbertInt_abs_max ...
        ... % Determinant
		... % умножители определителя 2x2 
		DetM_2x2_multiplier_total_abs_max, ...
		... % сумматоры определителя 2x2 
		Det2x2_sum_abs_max, ...
		... % умножители определителя 3х3
		Mult_DetM_3x3_array_max, ...
		... % разрядность умножителей определителя 3x3
		Mult_DetM_3x3_array_mult_total_width_max, ...
		... % пресумматоры определителя 3х3
		DetM_3x3_int_pre_sum_array_max, ...
		... % разрядность пресумматоров определителя 3x3
		DetM_3x3_int_pre_sum_width_total_max, ...
		... % сумматоры определителя 3х3
		DetM_3x3_int_sum_array_max, ...
		... % разрядность сумматоров определителя 3x3
		DetM_3x3_int_sum_width_total_max, ...
		... % умножители определителя 4х4
		DetM_4x4_int_mult_array_max, ...
		... % разрядность умножителей определителя 4х4
		DetM_4x4_int_mult_width_total_max, ...
		... % пресумматоры определителя 4х4
		DetM_4x4_int_pre_sum_array_max, ...
		... % разрядность пресумматоров определителя 4х4
		DetM_4x4_int_pre_sum_width_total_max, ...
		... % сумматоры определителя 4х4
		DetM_4x4_int_sum_array_max, ...
		... % разрядность сумматоров определителя 4х4
		DetM_4x4_int_sum_array_width_total_max, ...
		... % умножители определителя 5х5
		DetM_5x5_int_mult_array_max, ...
		... % разрядность умножителей 5x5
		DetM_5x5_int_mult_array_width_total_max, ...
		... % пресумматоры1 определителя 5х5
		DetM_5x5_int_pre_sum1_array_max, ...
		... % разрядность пресумматоров1 определителя 5х5
		DetM_5x5_int_pre_sum1_array_width_total_max, ...
		... % пресумматор2 определителя 5х5
		DetM_5x5_int_sum3_abs_max, ...
		... % разрядность пресумматора2 определителя 5х5
		DetM_5x5_int_sum3_width_total_max, ...
		... % сумматор определителя 5х5
		DetM_5x5_int_abs_max, ...
		... % разрядность сумматора определителя 5х5
		DetM_5x5_int_width_total_max, ...
        ... % начальный определитель
        Det_x3_int_max, ...
        ... % выход делителя
        Divide_max, ...
        ... % умножители адаптивного фильтра
        Adaptive_filter_mult_array_max, ...
        ... % разрядность умножителей адаптивного фильтра
        Adaptive_filter_mult_total_width, ...
        ... % сумматоры адаптивного фильтра
        Adaptive_filter_sum_array_max, ....
        ... % разрядность сумматоров адаптивного фильтра
        Adaptive_filter_sum_total_width ...
	] = adc_calibration(sim_options, adc_input, s_to_subadc, s_after_subadc);


    % Записываем значения каждого фильтра
    for j = 1:sim_options.N
            % выбираем максимальное значение сигнала умножителей 
            % полосового фильтра АЦП0
            if (golden_mult_max(j,i) < golden_mult(j,i))
                golden_mult_max(j,i) = golden_mult(j,i); 
            end
            % выбираем максимальное значение разрядности умножителя
            % полосового фильтра АЦП0
            if (golden_total_width_mult_max(j,i) < golden_width_total_mult(j,i))
                golden_total_width_mult_max(j,i) = golden_width_total_mult(j,i); 
            end
        end


    for i = 1:sim_options.M-1

        if (y_fractional_outInt_abs_max_in_cycle(i,1) < y_fractional_outInt_abs_max(i,1))
            y_fractional_outInt_abs_max_in_cycle(i,1) = y_fractional_outInt_abs_max(i,1); 
        end

        if (ymi_HilbertInt_abs_max_in_cycle(i,1) < ymi_HilbertInt_abs_max(i))
            ymi_HilbertInt_abs_max_in_cycle(i,1) = ymi_HilbertInt_abs_max(i); 
        end

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
        %% Determinante 2x2
        for k = 1:sim_options.num_det2x2*2
            % Поиск максимального значения умножителей в определителе 2х2
            % 2х2
            if (DetM_2x2_multiplier_total_abs_max_in_cycle(k,i) < DetM_2x2_multiplier_total_abs_max(k,i))
                DetM_2x2_multiplier_total_abs_max_in_cycle(k,i) = DetM_2x2_multiplier_total_abs_max(k,i);
            end
        end
		for k = 1:sim_options.num_det2x2
            % Поиск максимального значения сумматоров в определителе 2х2
            % 2х2
            if (Det2x2_sum_abs_max_in_cycle(k,i) < Det2x2_sum_abs_max(k,i))
                Det2x2_sum_abs_max_in_cycle(k,i) = Det2x2_sum_abs_max(k,i); 
            end
        end
		%% Определитель 3х3 умножители
        for k = 1:sim_options.num_det2x2*3
            if (Mult_DetM_3x3_array_max_in_cycle(k,i) < Mult_DetM_3x3_array_max(k,i))
                Mult_DetM_3x3_array_max_in_cycle(k,i) = Mult_DetM_3x3_array_max(k,i); 
            end
            if (Mult_DetM_3x3_array_mult_total_width_max_in_cycle(k,i) < Mult_DetM_3x3_array_mult_total_width_max(k,i))
                Mult_DetM_3x3_array_mult_total_width_max_in_cycle(k,i) = Mult_DetM_3x3_array_mult_total_width_max(k,i); 
            end
        end
		for k = 1:sim_options.num_det2x2
            if (DetM_3x3_int_pre_sum_array_max_in_cycle(k,i) < DetM_3x3_int_pre_sum_array_max(k,i))
                DetM_3x3_int_pre_sum_array_max_in_cycle(k,i) = DetM_3x3_int_pre_sum_array_max(k,i); 
            end
            if (DetM_3x3_int_pre_sum_width_total_max_in_cycle(k,i) < DetM_3x3_int_pre_sum_array_max(k,i))
                DetM_3x3_int_pre_sum_width_total_max_in_cycle(k,i) = DetM_3x3_int_pre_sum_array_max(k,i); 
            end
            if (DetM_3x3_int_pre_sum_width_total_max_in_cycle(k,i) < DetM_3x3_int_pre_sum_width_total_max(k,i))
                DetM_3x3_int_pre_sum_width_total_max_in_cycle(k,i) = DetM_3x3_int_pre_sum_width_total_max(k,i); 
            end
            if (DetM_3x3_int_sum_array_max_in_cycle(k,i) < DetM_3x3_int_sum_array_max(k,i))
                DetM_3x3_int_sum_array_max_in_cycle(k,i) = DetM_3x3_int_sum_array_max(k,i); 
            end
            if (DetM_3x3_int_sum_width_total_max_in_cycle(k,i) < DetM_3x3_int_sum_width_total_max(k,i))
                DetM_3x3_int_sum_width_total_max_in_cycle(k,i) = DetM_3x3_int_sum_width_total_max(k,i); 
            end
        end
		for k = 1:sim_options.num_det2x2*2
            if (DetM_4x4_int_mult_array_max_in_cycle(k,i) < DetM_4x4_int_mult_array_max(k,i))
                DetM_4x4_int_mult_array_max_in_cycle(k,i) = DetM_4x4_int_mult_array_max(k,i); 
            end
            if (DetM_4x4_int_mult_width_total_max_in_cycle(k,i) < DetM_4x4_int_mult_width_total_max(k,i))
                DetM_4x4_int_mult_width_total_max_in_cycle(k,i) = DetM_4x4_int_mult_width_total_max(k,i); 
            end
        end
		for k = 1:sim_options.num_det2x2
            if (DetM_4x4_int_pre_sum_array_max_in_cycle(k,i) < DetM_4x4_int_pre_sum_array_max(k,i))
                DetM_4x4_int_pre_sum_array_max_in_cycle(k,i) = DetM_4x4_int_pre_sum_array_max(k,i); 
            end
        end
		for k = 1:sim_options.Size_matrix
            if (DetM_4x4_int_sum_array_max_in_cycle(k,i) < DetM_4x4_int_sum_array_max(k,i))
                DetM_4x4_int_sum_array_max_in_cycle(k,i) = DetM_4x4_int_sum_array_max(k,i); 
            end
            if (DetM_5x5_int_mult_array_max_in_cycle(k,i) < DetM_5x5_int_mult_array_max(k,i))
                DetM_5x5_int_mult_array_max_in_cycle(k,i) = DetM_5x5_int_mult_array_max(k,i); 
            end
        end
		for k = 1:2
            if (DetM_5x5_int_pre_sum1_array_max_in_cycle(k,i) < DetM_5x5_int_pre_sum1_array_max(k,i))
                DetM_5x5_int_pre_sum1_array_max_in_cycle(k,i) = DetM_5x5_int_pre_sum1_array_max(k,i); 
            end
        end
        if (DetM_5x5_int_sum3_abs_max_in_cycle(1,i) < DetM_5x5_int_sum3_abs_max(1,i))
            DetM_5x5_int_sum3_abs_max_in_cycle(1,i) = DetM_5x5_int_sum3_abs_max(1,i); 
        end
		if (DetM_5x5_int_abs_max_in_cycle(1,i) < DetM_5x5_int_abs_max(i))
            DetM_5x5_int_abs_max_in_cycle(1,i) = DetM_5x5_int_abs_max(i); 
        end

		if (Det_x3_int_max_in_cycle(1,i) < Det_x3_int_max(i))
            Det_x3_int_max_in_cycle(1,i) = Det_x3_int_max(i); 
        end
        %% Делитель
        % выход делителя
		if (Divide_max_in_cycle(1,i) < Divide_max(i))
            Divide_max_in_cycle(1,i) = Divide_max(i); 
        end
        %% Адаптивный фильтр
		for k = 1:sim_options.Size_matrix
            if (Adaptive_filter_mult_array_max_in_cycle(k,i) < Adaptive_filter_mult_array_max(k,i))
                Adaptive_filter_mult_array_max_in_cycle(k,i) = Adaptive_filter_mult_array_max(k,i); 
            end
            if (Adaptive_filter_mult_total_width_in_cycle(k,i) < Adaptive_filter_mult_total_width(k,i))
                Adaptive_filter_mult_total_width_in_cycle(k,i) = Adaptive_filter_mult_total_width(k,i); 
            end
            if (Adaptive_filter_sum_array_max_in_cycle(k,i) < Adaptive_filter_sum_array_max(k,i))
                Adaptive_filter_sum_array_max_in_cycle(k,i) = Adaptive_filter_sum_array_max(k,i); 
            end
            if (Adaptive_filter_sum_total_width_in_cycle(k,i) < Adaptive_filter_sum_total_width(k,i))
                Adaptive_filter_sum_total_width_in_cycle(k,i) = Adaptive_filter_sum_total_width(k,i); 
            end
        end
    end

    %% SFDR
    figure(8);
    subplot(5,1,1);
    sfdr(s_to_subadc(1:length(x_after_adc)), sim_options.Fs/sim_options.Inter);
    subplot(5,1,2);
    sfdr(s_after_subadc(1:length(x_after_adc)), sim_options.Fs/sim_options.Inter);
    subplot(5,1,3);
    sfdr(x_after_adc(1:length(x_after_adc)), sim_options.Fs/sim_options.Inter);
    subplot(5,1,4);
    sfdr(x_after_adc_double(1:length(x_after_adc_double)), sim_options.Fs/sim_options.Inter);
    subplot(5,1,5);
    sfdr(x_after_adc_int(1:length(x_after_adc)), sim_options.Fs/sim_options.Inter);
    %% SNR
    figure(9);
    subplot(5,1,1);
    snr(s_to_subadc(1:length(x_after_adc)), sim_options.Fs/sim_options.Inter);
    subplot(5,1,2);
    snr(s_after_subadc(1:length(x_after_adc)), sim_options.Fs/sim_options.Inter);
    subplot(5,1,3);
    snr(x_after_adc(1:length(x_after_adc)), sim_options.Fs/sim_options.Inter);
    subplot(5,1,4);
    snr(x_after_adc_double(100:length(x_after_adc)), sim_options.Fs/sim_options.Inter);
    subplot(5,1,5);
    snr(x_after_adc_int(100:length(x_after_adc)), sim_options.Fs/sim_options.Inter);
    %% 
    snr_in_int(num) = snr(double(s_after_subadc), sim_options.Fs/sim_options.Inter);
    snr_output_lu(num) = snr(x_after_adc, sim_options.Fs/sim_options.Inter);
    snr_output_double(num) = snr(x_after_adc_double, sim_options.Fs/sim_options.Inter);
    snr_output_int(num) = snr(x_after_adc_int, sim_options.Fs/sim_options.Inter);

    sfdr_in_int(num) = sfdr(double(s_after_subadc), sim_options.Fs/sim_options.Inter);
    sfdr_output_lu(num) = sfdr(x_after_adc, sim_options.Fs/sim_options.Inter);
    sfdr_output_double(num) = sfdr(x_after_adc_double, sim_options.Fs/sim_options.Inter);
    sfdr_output_int(num) = sfdr(x_after_adc_int, sim_options.Fs/sim_options.Inter);

    norm_freq(num) = sim_options.freq/(sim_options.Fs/sim_options.Inter/sim_options.M);
    num_array(:,num) = num;

    % freq
    sim_options.freq = sim_options.freq + sim_options.step; % frequency of fundamental tone
    % SNR
    sim_options.SNR = sim_options.SNR + sim_options.Step_of_SNR;
end


if sim_options.enable_mask == false

    %% Запись данных для полосового фильтра АЦП0
    % формируем таблицу максимальных значений полосового фильтра АЦП0
    T1 = table(golden_mult, golden_sum, 'RowNames', num_filter);
    writetable(T1,['src/width_txt/Максимальные_значения_полосового_фильтра_АЦП_0.xlsx'],'WriteRowNames',true);  

    % формируем таблицу максимальных разрядностей полосового фильтра АЦП0
    T2 = table(golden_total_width_mult_max, golden_total_width_sum_max, 'RowNames', num_filter);
    writetable(T2,['src/width_txt/Разрядность_элементов_полосового_фильтра_АЦП_0.xlsx'],'WriteRowNames',true); 

    % формируем таблицу максимальных значений выхода фильтра дробной задержки
    y_fractional_outInt = y_fractional_outInt_abs_max_in_cycle(:,1);
    T3 = table(y_fractional_outInt, 'RowNames', num_filter_out);
        writetable(T3,['src/width_txt/Максимальные_выходные_значения_фильтра_дробной_задержки_АЦП_1_3.xlsx'],'WriteRowNames',true); 

    % формируем таблицу максимальных значений выхода фильтра дробной задержки
    y_Hilbert_outInt = ymi_HilbertInt_abs_max_in_cycle(:,1);
    T4 = table(y_Hilbert_outInt, 'RowNames', num_filter_out);
       writetable(T4,['src/width_txt/Максимальные_выходные_значения_фильтра_Гилберта_АЦП_1_3.xlsx'],'WriteRowNames',true); 

    for i = 1:sim_options.M-1
        for j = 1:length(fractional_mult_max(:,i))
            width_mult(j,i) = define_of_width_int(fractional_mult_max(j,i), sim_options.int_size, sim_options.width_fractional);
            width_mult_h(j,i) = define_of_width_int(hilbert_mult_max(j,i), sim_options.int_size, sim_options.width_hilbert);
        end
        for j = 1:length(fractional_sum_max(:,i))-1
            width_sum(j,i) = define_of_width_int(fractional_sum_max(j,i), sim_options.int_size, sim_options.width_fractional);
            width_sum_h(j,i) = define_of_width_int(hilbert_sum_max(j,i), sim_options.int_size, sim_options.width_hilbert);
        end

        %% Запись данных для фильтра дробной задержки
        % запись макс. значений сигнала
        Multipliers_fractional = fractional_mult_max(:,i);
        Adders_fractional = fractional_sum_max(:,i);
        % формируем таблицу максимальных значений фильтра дробной задержки
        T1 = table(Multipliers_fractional, Adders_fractional, 'RowNames', num_filter);
        writetable(T1,['src/width_txt/Максимальные_значения_фильтра_дробной_задержки_АЦП_№' num2str(i) '.xlsx'],'WriteRowNames',true);  

        Multipliers_fractional_width = fractional_total_width_mult_max(:,i);
        Adders_fractional_width = fractional_total_width_sum_max(:,i);
       
        % формируем таблицу максимальных разрядностей фильтра дробной задержки
        T3 = table(Multipliers_fractional_width, Adders_fractional_width, 'RowNames', num_filter);
        writetable(T3,['src/width_txt/Разрядность_элементов_фильтра_дробной_задержки_АЦП_№' num2str(i) '.xlsx'],'WriteRowNames',true); 
        
        %% Запись данных для фильтра Гилберта
        % запись макс. значений сигнала
        Multipliers_hilbert = hilbert_mult_max(:,i);
        Adders_hilbert= hilbert_sum_max(:,i);
        % формируем таблицу фильтра дробной задержки
        T4 = table(Multipliers_hilbert, Adders_hilbert, 'RowNames', num_filter);
        writetable(T4,['src/width_txt/Максимальные_значения_фильтра_Гилберта_АЦП_№' num2str(i) '.xlsx'],'WriteRowNames',true);  

        Multipliers_hilbert_width = hilbert_total_width_mult_max(:,i);
        Adders_hilbert_width = hilbert_total_width_sum_max(:,i);
 
        % формируем таблицу максимальных разрядностей фильтра дробной задержки
        T5 = table(Multipliers_hilbert_width, Adders_hilbert_width, 'RowNames', num_filter);
        writetable(T5,['src/width_txt/Разрядность_элементов_фильтра_Гилберта_АЦП_№' num2str(i) '.xlsx'],'WriteRowNames',true); 

        %% Запись данных для определителя
        Multipliers_2x2 = DetM_2x2_multiplier_total_abs_max_in_cycle(:,i);
        Adders_2x2 = Det2x2_sum_abs_max_in_cycle(:,i);
        Multipliers_3x3 = Mult_DetM_3x3_array_max_in_cycle(:,i);
        Pre_sum_3x3 = DetM_3x3_int_pre_sum_array_max_in_cycle(:,i);
        Sum_3x3 = DetM_3x3_int_sum_array_max_in_cycle(:,i);
		Multipliers_4x4 = DetM_4x4_int_mult_array_max_in_cycle(:,i);
		Pre_sum_4x4 = DetM_4x4_int_pre_sum_array_max_in_cycle(:,i);
		Sum4x4 = DetM_4x4_int_sum_array_max_in_cycle(:,i);
		Mult5x5 = DetM_5x5_int_mult_array_max_in_cycle(:,i);
		Pre_sum1_5x5 = DetM_5x5_int_pre_sum1_array_max_in_cycle(:,i);
		Pre_sum2_5x5 = DetM_5x5_int_sum3_abs_max_in_cycle(:,i);
		Det_5x5 = DetM_5x5_int_abs_max_in_cycle(:,i);
        % формируем таблицу определителя
        T6 = table(Multipliers_2x2, Adders_2x2, Multipliers_3x3, Pre_sum_3x3, Sum_3x3, Multipliers_4x4, Pre_sum_4x4, Sum4x4, Mult5x5, Pre_sum1_5x5, Pre_sum2_5x5, Det_5x5, 'RowNames', num_determinante);
        writetable(T6,['src/width_txt/Максимальные_значения_определителя_АЦП_№' num2str(i) '.xlsx'],'WriteRowNames',true)  
        %% Запись данных для делителя
        Dividend = Det_5x5(1,1);
        Divisor = Det_x3_int_max_in_cycle(:,i);
        Quotient = Divide_max_in_cycle(:,i);
        % формируем таблицу делителя
        T7 = table(Dividend, Divisor, Quotient, 'RowNames', num_divide);
        writetable(T7,['src/width_txt/Максимальные_значения_делителя_АЦП_№' num2str(i) '.xlsx'],'WriteRowNames',true)  
        
        %% Запись данных для адаптивного фильтра
        Multipliers_adaptive = Adaptive_filter_mult_array_max_in_cycle(:,i);
        Sum_adaptive = Adaptive_filter_sum_array_max_in_cycle(:,i);
        % формируем таблицу адаптивного фильтра
        T8 = table(Multipliers_adaptive, Sum_adaptive, 'RowNames', num_adaptive);
        writetable(T8,['src/width_txt/Максимальные_значения_адаптивного_фильтра_АЦП_№' num2str(i) '.xlsx'],'WriteRowNames',true) 

        Multipliers_adaptive_total_width = Adaptive_filter_mult_total_width_in_cycle(:,i);
        Sum_adaptive_total_width = Adaptive_filter_sum_total_width_in_cycle(:,i);
        % формируем таблицу адаптивного фильтра
        T9 = table(Multipliers_adaptive_total_width, Sum_adaptive_total_width, 'RowNames', num_adaptive);
        writetable(T9,['src/width_txt/Разрядность_элементов_адаптивного_фильтра_АЦП_№' num2str(i) '.xlsx'],'WriteRowNames',true) 

    end
end

%% Итоговый график SNR и SFDR каждой итерации алгоритма
figure(10);
subplot(2,1,1)
plot(norm_freq, snr_in_int, '-o', norm_freq, snr_output_lu, '-o', norm_freq, snr_output_double, '-o', norm_freq, snr_output_int, '-o');
title('SNR')
xlabel('Нормированная частота') 
ylabel('SNR (dB)') 
legend({'Входной сигнал с ошибками int', 'Выходной сигнал матлаб функции LU', 'Выходной сигнал double', 'Выходной сигнал int'}, 'Location','northwest');
% 
subplot(2,1,2)
plot(norm_freq, sfdr_in_int, '-o', norm_freq, sfdr_output_lu, '-o', norm_freq, sfdr_output_double, '-o', norm_freq, sfdr_output_int, '-o');
title('SFDR (dB)')
xlabel('Нормированная частота') 
ylabel('SFDR (dB)') 
legend({'Входной сигнал с ошибками int', 'Выходной сигнал матлаб функции LU', 'Выходной сигнал double', 'Выходной сигнал int'}, 'Location','northwest');

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
% 
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