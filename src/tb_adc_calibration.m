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

num_filter = string((1:sim_options.N)');
num_determinante = string((1:sim_options.num_det2x2*3)');
num_adaptive = string((1:sim_options.Size_matrix*sim_options.Size_matrix)');
num_divide = string((1:1)');
%% Полосовой фильтр АЦП0
bandpass_mult_max = cast(zeros(sim_options.N,1), sim_options.int_size);
bandpass_sum_max = cast(zeros(sim_options.N,1), sim_options.int_size);
bandpass_total_width_mult_max = cast(zeros(sim_options.N,1), sim_options.int_size);
bandpass_total_width_sum_max = cast(zeros(sim_options.N,1), sim_options.int_size);
bandpass_width_mult = cast(zeros(sim_options.N,1), sim_options.int_size);
bandpass_width_sum = cast(zeros(sim_options.N,1), sim_options.int_size);
%% Фильтр дробной задержки
width_mult = cast(zeros(sim_options.N, sim_options.M-1), sim_options.int_size);
width_sum = cast(zeros(sim_options.N, sim_options.M-1), sim_options.int_size);
fractional_mult_max = cast(zeros(sim_options.N,sim_options.M-1), sim_options.int_size);
fractional_sum_max = cast(zeros(sim_options.N,sim_options.M-1), sim_options.int_size);
fractional_total_width_mult_max = cast(zeros(sim_options.N,sim_options.M-1), sim_options.int_size);
fractional_total_width_sum_max = cast(zeros(sim_options.N,sim_options.M-1), sim_options.int_size);
y_fractional_outInt_abs_max_in_cycle = cast(zeros(sim_options.N, sim_options.M-1), sim_options.int_size);
%% Фильтр Гилберта
width_mult_h = cast(zeros(sim_options.N, sim_options.M-1), sim_options.int_size);
width_sum_h = cast(zeros(sim_options.N, sim_options.M-1), sim_options.int_size);
hilbert_mult_max = cast(zeros(sim_options.N,sim_options.M-1), sim_options.int_size);
hilbert_sum_max = cast(zeros(sim_options.N,sim_options.M-1), sim_options.int_size);
hilbert_total_width_mult_max = cast(zeros(sim_options.N,sim_options.M-1), sim_options.int_size);
hilbert_total_width_sum_max = cast(zeros(sim_options.N,sim_options.M-1), sim_options.int_size);
ymi_HilbertInt_abs_max_in_cycle = cast(zeros(sim_options.N,sim_options.M-1), sim_options.int_size);
%%
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


Multipliers_2x2 = cast(zeros(sim_options.num_det2x2*3,1), sim_options.type_2x2_det);
Adders_2x2 = cast(zeros(sim_options.num_det2x2*3, 1), sim_options.type_2x2_det);
Multipliers_3x3 = cast(zeros(sim_options.num_det2x2*3, 1), sim_options.type_3x3_det);
Pre_sum_3x3 = cast(zeros(sim_options.num_det2x2*3, 1), sim_options.type_3x3_det);
Sum_3x3 = cast(zeros(sim_options.num_det2x2*3, 1), sim_options.type_3x3_det);

Multipliers_4x4 = cast(zeros(sim_options.num_det2x2*3, 1), sim_options.type_4x4_det);
Pre_sum_4x4 = cast(zeros(sim_options.num_det2x2*3, 1), sim_options.type_4x4_det);
Sum4x4 = cast(zeros(sim_options.num_det2x2*3, 1), sim_options.type_4x4_det);

Mult5x5 = cast(zeros(sim_options.num_det2x2*3, 1), sim_options.type_5x5_det);
Pre_sum1_5x5 = cast(zeros(sim_options.num_det2x2*3, 1), sim_options.type_5x5_det);
Pre_sum2_5x5 = cast(zeros(sim_options.num_det2x2*3, 1), sim_options.type_5x5_det);
Det_5x5 = cast(zeros(sim_options.num_det2x2*3, 1), sim_options.type_5x5_det);
%%
width_mult_det2x2 = cast(zeros(length(num_determinante),1), sim_options.type_2x2_det);
width_sum_det2x2 = cast(zeros(length(num_determinante),1), sim_options.type_2x2_det);

width_mult_det3x3 = cast(zeros(length(num_determinante),1), sim_options.type_3x3_det);
width_presum_det3x3 = cast(zeros(length(num_determinante),1), sim_options.type_3x3_det);
width_sum_det3x3 = cast(zeros(length(num_determinante),1), sim_options.type_3x3_det);

width_mult_det4x4 = cast(zeros(length(num_determinante),1), sim_options.type_4x4_det);
width_presum_det4x4 = cast(zeros(length(num_determinante),1), sim_options.type_4x4_det);
width_sum_det4x4 = cast(zeros(length(num_determinante),1), sim_options.type_4x4_det);

width_mult_det5x5 = cast(zeros(length(num_determinante),1), sim_options.type_5x5_det);
width_presum1_det5x5 = cast(zeros(length(num_determinante),1), sim_options.type_5x5_det);
width_presum2_det5x5 = cast(zeros(length(num_determinante),1), sim_options.type_5x5_det);
width_sum_det5x5 = cast(zeros(length(num_determinante),1), sim_options.type_5x5_det);

% начальный определитель
Det_x3_int_max_in_cycle = cast(zeros(1,sim_options.M-1), sim_options.type_5x5_det);
% выход делителя
Divide_max_in_cycle = cast(zeros(1,sim_options.M-1), sim_options.type_divide_out);

%% Адаптивный фильтр
Adaptive_filter_mult_array_max_in_cycle = cast(zeros(sim_options.Size_matrix*sim_options.Size_matrix,sim_options.M-1), sim_options.type_mult_in_adaptive_filter);
Adaptive_filter_mult_total_width_in_cycle = cast(zeros(sim_options.Size_matrix*sim_options.Size_matrix,sim_options.M-1), sim_options.type_mult_in_adaptive_filter);
Adaptive_filter_sum_array_max_in_cycle = cast(zeros(sim_options.Size_matrix*sim_options.Size_matrix, 3, sim_options.M-1), sim_options.type_add_in_adaptive_filter);
Adaptive_filter_sum_total_width_in_cycle = cast(zeros(sim_options.Size_matrix*sim_options.Size_matrix, 3, sim_options.M-1), sim_options.type_add_in_adaptive_filter);

Multipliers_adaptive = cast(zeros(sim_options.Size_matrix*sim_options.Size_matrix, 1), sim_options.type_mult_in_adaptive_filter);
Sum1_adaptive = cast(zeros(sim_options.Size_matrix*sim_options.Size_matrix, 1), sim_options.type_add_in_adaptive_filter);
Sum2_adaptive = cast(zeros(sim_options.Size_matrix*sim_options.Size_matrix, 1), sim_options.type_add_in_adaptive_filter);
Sum3_adaptive =cast(zeros(sim_options.Size_matrix*sim_options.Size_matrix, 1), sim_options.type_add_in_adaptive_filter);
		
width_mult_adaptive = cast(zeros(length(num_adaptive),1), sim_options.type_mult_in_adaptive_filter);
width_sum1_adaptive = cast(zeros(length(num_adaptive),1), sim_options.type_add_in_adaptive_filter);
width_sum2_adaptive = cast(zeros(length(num_adaptive),1), sim_options.type_add_in_adaptive_filter);
width_sum3_adaptive = cast(zeros(length(num_adaptive),1), sim_options.type_add_in_adaptive_filter);

            
for num = 1:sim_options.num_cycles

    % Функция генерации сигналов для АЦП
    [s_to_subadc, adc_input, s_after_subadc, sim_options.Z] = gen_oversampled_signal(sim_options);

    % Основная функция калибровки АЦП
    [x_after_adc, x_after_adc_double, x_after_adc_int, ...
        ... % Значения полосового фильтра ADC0, ...
        filter_max_width_out_golden, ...
        ... % Значения фильтров дробной задержки
        filter_max_width_out_fractional, ...
        y_fractional_outInt_abs_max, ...
        ... % Значения фильтров Гилберта
        filter_max_width_out_hilbert, ... % hilbert_mult, hilbert_sum, hilbert_width_total_mult, hilbert_width_total_sum, ymi_HilbertInt_abs_max, ...
        ymi_HilbertInt_abs_max, ...
        ... % Determinant
        determinate_struct, ...
        ... % выход делителя
        Divide_max, ...
        adaptive_filter_struct_max ...
	] = adc_calibration(sim_options, adc_input, s_to_subadc, s_after_subadc);


   % save (sprintf(num2str(clock) + ".mat"));
   % load ('2025             12              8             17             14         49.333.mat');

   
   %% Поиск максимума 

    % Полосовой фильтр АЦП0
    % Записываем значения полосового фильтра АЦП0
	% умножители
    for j = 1:sim_options.N
        % выбираем максимальное значение сигнала умножителей 
        % полосового фильтра АЦП0
        if (bandpass_mult_max(j) < filter_max_width_out_golden.mult_max(j))
            bandpass_mult_max(j) = filter_max_width_out_golden.mult_max(j); 
        end
		
        % выбираем максимальное значение разрядности умножителя
        % полосового фильтра АЦП0
        if (bandpass_total_width_mult_max(j) < filter_max_width_out_golden.width_total_mult_max(j))
            bandpass_total_width_mult_max(j) = filter_max_width_out_golden.width_total_mult_max(j); 
        end
    end
	
	% сумматоры
	for j = 1:sim_options.N-1
        % выбираем максимальное значение сигнала умножителей 
        % полосового фильтра АЦП0
        if (bandpass_sum_max(j) < filter_max_width_out_golden.sum_max(j))
            bandpass_sum_max(j) = filter_max_width_out_golden.sum_max(j); 
        end
		
        % выбираем максимальное значение разрядности сумматоров
        % полосового фильтра АЦП0
        if (bandpass_total_width_sum_max(j) < filter_max_width_out_golden.width_total_sum_max(j))
            bandpass_total_width_sum_max(j) = filter_max_width_out_golden.width_total_sum_max(j); 
        end
    end
	
	%%

    for i = 1:sim_options.M-1
	
        if (y_fractional_outInt_abs_max_in_cycle(1,i) < y_fractional_outInt_abs_max(i))
            y_fractional_outInt_abs_max_in_cycle(1,i) = y_fractional_outInt_abs_max(i); 
        end
        if (ymi_HilbertInt_abs_max_in_cycle(1,i) < ymi_HilbertInt_abs_max(i))
            ymi_HilbertInt_abs_max_in_cycle(1,i) = ymi_HilbertInt_abs_max(i); 
        end

		% умножители фильтров
        for j = 1:sim_options.N
            % выбираем максимальное значение сигнала умножителей фильтра
            % дробной задержки
            if (fractional_mult_max(j,i) < filter_max_width_out_fractional(i).mult_max(j))
                fractional_mult_max(j,i) = filter_max_width_out_fractional(i).mult_max(j); 
            end
            % выбираем максимальное значение разрядности умножителя фильтра
            % дробной задержки
            if (fractional_total_width_mult_max(j,i) < filter_max_width_out_fractional(i).width_total_mult_max(j))
                fractional_total_width_mult_max(j,i) = filter_max_width_out_fractional(i).width_total_mult_max(j); 
            end
            %%
            % выбираем максимальное значение сигнала умножителей фильтра
            % Гилберта
            if (hilbert_mult_max(j,i) < filter_max_width_out_hilbert(i).mult_max(j))
                hilbert_mult_max(j,i) = filter_max_width_out_hilbert(i).mult_max(j); 
            end
            % выбираем максимальное значение разрядности умножителя фильтра
            % Гилберта
            if (hilbert_total_width_mult_max(j,i) < filter_max_width_out_hilbert(i).width_total_mult_max(j))
                hilbert_total_width_mult_max(j,i) = filter_max_width_out_hilbert(i).width_total_mult_max(j); 
            end
        end
		% сумматоры фильтров
        for j = 1:sim_options.N-1
            % выбираем максимальное значение сигнала сумматоров
            % фильтра дробной задержки 
            if (fractional_sum_max(j,i) < filter_max_width_out_fractional(i).sum_max(j))
                fractional_sum_max(j,i) = filter_max_width_out_fractional(i).sum_max(j); 
            end
            % выбираем максимальное значение разрядности сумматоров
            % фильтра дробной задержки
            if (fractional_total_width_sum_max(j,i) < filter_max_width_out_fractional(i).width_total_sum_max(j))
                fractional_total_width_sum_max(j,i) = filter_max_width_out_fractional(i).width_total_sum_max(j); 
            end
            %%
            % выбираем максимальное значение сигнала сумматоров
            % фильтра Гилберта
            if (hilbert_sum_max(j,i) < filter_max_width_out_hilbert(i).sum_max(j))
                hilbert_sum_max(j,i) = filter_max_width_out_hilbert(i).sum_max(j); 
            end
            % выбираем максимальное значение разрядности сумматоров
            % фильтра Гилберта
            if (hilbert_total_width_sum_max(j,i) < filter_max_width_out_hilbert(i).width_total_sum_max(j))
                hilbert_total_width_sum_max(j,i) = filter_max_width_out_hilbert(i).width_total_sum_max(j); 
            end
        end


        %% Determinante 2x2
        for k = 1:sim_options.num_det2x2*2
            % Поиск максимального значения умножителей в определителе 2х2
            % 2х2
            if (DetM_2x2_multiplier_total_abs_max_in_cycle(k,i) < determinate_struct(i).DetM_2x2_multiplier_total_abs_max(k))
                DetM_2x2_multiplier_total_abs_max_in_cycle(k,i) = determinate_struct(i).DetM_2x2_multiplier_total_abs_max(k);
            end
        end
		for k = 1:sim_options.num_det2x2
            % Поиск максимального значения сумматоров в определителе 2х2
            % 2х2
            if (Det2x2_sum_abs_max_in_cycle(k,i) < determinate_struct(i).Det2x2_sum_abs_max(k))
                Det2x2_sum_abs_max_in_cycle(k,i) = determinate_struct(i).Det2x2_sum_abs_max(k); 
            end
        end
		%% Определитель 3х3 умножители
        for k = 1:sim_options.num_det2x2*3
            if (Mult_DetM_3x3_array_max_in_cycle(k,i) < determinate_struct(i).Mult_DetM_3x3_array_max(k))
                Mult_DetM_3x3_array_max_in_cycle(k,i) = determinate_struct(i).Mult_DetM_3x3_array_max(k); 
            end
            if (Mult_DetM_3x3_array_mult_total_width_max_in_cycle(k,i) < determinate_struct(i).Mult_DetM_3x3_array_mult_total_width_max(k))
                Mult_DetM_3x3_array_mult_total_width_max_in_cycle(k,i) = determinate_struct(i).Mult_DetM_3x3_array_mult_total_width_max(k); 
            end
        end
		for k = 1:sim_options.num_det2x2
            if (DetM_3x3_int_pre_sum_array_max_in_cycle(k,i) < determinate_struct(i).DetM_3x3_int_pre_sum_array_max(k))
                DetM_3x3_int_pre_sum_array_max_in_cycle(k,i) = determinate_struct(i).DetM_3x3_int_pre_sum_array_max(k); 
            end
            if (DetM_3x3_int_pre_sum_width_total_max_in_cycle(k,i) < determinate_struct(i).DetM_3x3_int_pre_sum_array_max(k))
                DetM_3x3_int_pre_sum_width_total_max_in_cycle(k,i) = determinate_struct(i).DetM_3x3_int_pre_sum_array_max(k); 
            end
            if (DetM_3x3_int_pre_sum_width_total_max_in_cycle(k,i) < determinate_struct(i).DetM_3x3_int_pre_sum_width_total_max(k))
                DetM_3x3_int_pre_sum_width_total_max_in_cycle(k,i) = determinate_struct(i).DetM_3x3_int_pre_sum_width_total_max(k); 
            end
            if (DetM_3x3_int_sum_array_max_in_cycle(k,i) < determinate_struct(i).DetM_3x3_int_sum_array_max(k))
                DetM_3x3_int_sum_array_max_in_cycle(k,i) = determinate_struct(i).DetM_3x3_int_sum_array_max(k); 
            end
            if (DetM_3x3_int_sum_width_total_max_in_cycle(k,i) < determinate_struct(i).DetM_3x3_int_sum_width_total_max(k))
                DetM_3x3_int_sum_width_total_max_in_cycle(k,i) = determinate_struct(i).DetM_3x3_int_sum_width_total_max(k); 
            end
        end
		for k = 1:sim_options.num_det2x2*2
            if (DetM_4x4_int_mult_array_max_in_cycle(k,i) < determinate_struct(i).DetM_4x4_int_mult_array_max(k))
                DetM_4x4_int_mult_array_max_in_cycle(k,i) = determinate_struct(i).DetM_4x4_int_mult_array_max(k); 
            end
            if (DetM_4x4_int_mult_width_total_max_in_cycle(k,i) < determinate_struct(i).DetM_4x4_int_mult_width_total_max(k))
                DetM_4x4_int_mult_width_total_max_in_cycle(k,i) = determinate_struct(i).DetM_4x4_int_mult_width_total_max(k); 
            end
        end
		for k = 1:sim_options.num_det2x2
            if (DetM_4x4_int_pre_sum_array_max_in_cycle(k,i) < determinate_struct(i).DetM_4x4_int_pre_sum_array_max(k))
                DetM_4x4_int_pre_sum_array_max_in_cycle(k,i) = determinate_struct(i).DetM_4x4_int_pre_sum_array_max(k); 
            end
        end
		for k = 1:sim_options.Size_matrix
            if (DetM_4x4_int_sum_array_max_in_cycle(k,i) < determinate_struct(i).DetM_4x4_int_sum_array_max(k))
                DetM_4x4_int_sum_array_max_in_cycle(k,i) = determinate_struct(i).DetM_4x4_int_sum_array_max(k); 
            end
            if (DetM_5x5_int_mult_array_max_in_cycle(k,i) < determinate_struct(i).DetM_5x5_int_mult_array_max(k))
                DetM_5x5_int_mult_array_max_in_cycle(k,i) = determinate_struct(i).DetM_5x5_int_mult_array_max(k); 
            end
        end
		for k = 1:2
            if (DetM_5x5_int_pre_sum1_array_max_in_cycle(k,i) < determinate_struct(i).DetM_5x5_int_pre_sum1_array_max(k))
                DetM_5x5_int_pre_sum1_array_max_in_cycle(k,i) = determinate_struct(i).DetM_5x5_int_pre_sum1_array_max(k); 
            end
        end
        if (DetM_5x5_int_sum3_abs_max_in_cycle(1,i) < determinate_struct(i).DetM_5x5_int_sum3_abs_max)
            DetM_5x5_int_sum3_abs_max_in_cycle(1,i) = determinate_struct(i).DetM_5x5_int_sum3_abs_max; 
        end
		if (DetM_5x5_int_abs_max_in_cycle(1,i) < determinate_struct(i).DetM_5x5_int_abs_max)
            DetM_5x5_int_abs_max_in_cycle(1,i) = determinate_struct(i).DetM_5x5_int_abs_max; 
        end

		if (Det_x3_int_max_in_cycle(1,i) < determinate_struct(i).Det_x3_int_max)
            Det_x3_int_max_in_cycle(1,i) = determinate_struct(i).Det_x3_int_max; 
        end
        %% Делитель
        % выход делителя
		if (Divide_max_in_cycle(1,i) < Divide_max(i))
            Divide_max_in_cycle(1,i) = Divide_max(i); 
        end
        %% Адаптивный фильтр
        % умножителм
        for k = 1:sim_options.Size_matrix*sim_options.Size_matrix
            if (Adaptive_filter_mult_array_max_in_cycle(k,i) < adaptive_filter_struct_max(i).Adaptive_filter_mult_array_max(k))
                Adaptive_filter_mult_array_max_in_cycle(k,i) = adaptive_filter_struct_max(i).Adaptive_filter_mult_array_max(k); 
            end
            if (Adaptive_filter_mult_total_width_in_cycle(k,i) < adaptive_filter_struct_max(i).Adaptive_filter_mult_total_width(k))
                Adaptive_filter_mult_total_width_in_cycle(k,i) = adaptive_filter_struct_max(i).Adaptive_filter_mult_total_width(k); 
            end
        end
        % сумматоры
        for t = 1:3
		    for k = 1:sim_options.Size_matrix*2
                if (Adaptive_filter_sum_array_max_in_cycle(k,t,i) < adaptive_filter_struct_max(i).Adaptive_filter_sum_array_max(k,t))
                    Adaptive_filter_sum_array_max_in_cycle(k,t,i) = adaptive_filter_struct_max(i).Adaptive_filter_sum_array_max(k,t); 
                end
                if (Adaptive_filter_sum_total_width_in_cycle(k,t,i) < adaptive_filter_struct_max(i).Adaptive_filter_sum_total_width(k,t))
                    Adaptive_filter_sum_total_width_in_cycle(k,t,i) = adaptive_filter_struct_max(i).Adaptive_filter_sum_total_width(k,t); 
                end
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
    T1 = table(bandpass_mult_max, bandpass_sum_max, 'RowNames', num_filter);
    writetable(T1,['src/width_txt/Максимальные_значения_полосового_фильтра_АЦП_0.xlsx'],'WriteRowNames',true); 
	
	% формируем таблицу максимальных разрядностей полосового фильтра АЦП0
    T2 = table(bandpass_total_width_mult_max, bandpass_total_width_sum_max, 'RowNames', num_filter);
    writetable(T2,['src/width_txt/Разрядность_элементов_полосового_фильтра_АЦП_0.xlsx'],'WriteRowNames',true); 
	
	% Запись разрядностей макс.значений
	for i = 1:sim_options.N
		bandpass_width_mult(i) = define_of_width_int(bandpass_mult_max(i), sim_options.int_size, sim_options.width_fractional);
	end
	for i = 1:sim_options.N-1
		bandpass_width_sum(i) = define_of_width_int(bandpass_sum_max(i), sim_options.int_size, sim_options.width_fractional);
	end
	
	% формируем таблицу максимальных разрядностей полосового фильтра АЦП0
    T3 = table(bandpass_width_mult, bandpass_width_sum, 'RowNames', num_filter);
    writetable(T3,['src/width_txt/Разрядность_максимальных_значений_эталонного_фильтра.xlsx'],'WriteRowNames',true);
	

    for i = 1:sim_options.M-1
        for j = 1:length(fractional_mult_max(:,i))
            width_mult(j,i) = define_of_width_int(fractional_mult_max(j,i), sim_options.int_size, sim_options.width_fractional);
            width_mult_h(j,i) = define_of_width_int(hilbert_mult_max(j,i), sim_options.int_size, sim_options.width_hilbert);
        end
        for j = 1:length(fractional_sum_max(:,i))-1
            width_sum(j,i) = define_of_width_int(fractional_sum_max(j,i), sim_options.int_size, sim_options.width_fractional);
            width_sum_h(j,i) = define_of_width_int(hilbert_sum_max(j,i), sim_options.int_size, sim_options.width_hilbert);
        end

        mult_fractional = width_mult(:,i);
        sum_fractional = width_sum(:,i);

        mult_hilbert = width_mult_h(:,i);
        sum_hilbert = width_sum_h(:,i);

        %% Запись данных для фильтра дробной задержки
        % запись макс. значений сигнала
        Multipliers_fractional = fractional_mult_max(:,i);
        Adders_fractional = fractional_sum_max(:,i);
		y_fractional_outInt = y_fractional_outInt_abs_max_in_cycle(:,i);
        % формируем таблицу максимальных значений фильтра дробной задержки
        T4 = table(Multipliers_fractional, Adders_fractional, y_fractional_outInt, 'RowNames', num_filter);
        writetable(T4,['src/width_txt/Максимальные_значения_фильтра_дробной_задержки_АЦП_№' num2str(i) '.xlsx'],'WriteRowNames',true);  

        Multipliers_fractional_width = fractional_total_width_mult_max(:,i);
        Adders_fractional_width = fractional_total_width_sum_max(:,i);
        % формируем таблицу максимальных разрядностей фильтра дробной задержки
        T5 = table(Multipliers_fractional_width, Adders_fractional_width, 'RowNames', num_filter);
        writetable(T5,['src/width_txt/Разрядность_элементов_фильтра_дробной_задержки_АЦП_№' num2str(i) '.xlsx'],'WriteRowNames',true); 
		
        % формируем таблицу разрядностей макс.значений фильтра дробной задержки
        T6 = table(mult_fractional, sum_fractional, 'RowNames', num_filter);
        writetable(T6,['src/width_txt/Разрядность_максимальных_значений_фильтра_дробной_задержки_АЦП_№' num2str(i) '.xlsx'],'WriteRowNames',true);
		
   
        %% Запись данных для фильтра Гилберта
        % запись макс. значений сигнала
        Multipliers_hilbert = hilbert_mult_max(:,i);
        Adders_hilbert = hilbert_sum_max(:,i);
		y_Hilbert_outInt = ymi_HilbertInt_abs_max_in_cycle(:,i);
        % формируем таблицу фильтра дробной задержки
        T7 = table(Multipliers_hilbert, Adders_hilbert, y_Hilbert_outInt, 'RowNames', num_filter);
        writetable(T7,['src/width_txt/Максимальные_значения_фильтра_Гилберта_АЦП_№' num2str(i) '.xlsx'],'WriteRowNames',true);  

        Multipliers_hilbert_width = hilbert_total_width_mult_max(:,i);
        Adders_hilbert_width = hilbert_total_width_sum_max(:,i);
        % формируем таблицу максимальных разрядностей фильтра дробной задержки
        T8 = table(Multipliers_hilbert_width, Adders_hilbert_width, 'RowNames', num_filter);
        writetable(T8,['src/width_txt/Разрядность_элементов_фильтра_Гилберта_АЦП_№' num2str(i) '.xlsx'],'WriteRowNames',true); 

		% формируем таблицу разрядностей макс.значений фильтра Гилберта
        T9 = table(mult_hilbert, sum_hilbert, 'RowNames', num_filter);
        writetable(T9,['src/width_txt/Разрядность_максимальных_значений_фильтра_Гилберта_АЦП_№' num2str(i) '.xlsx'],'WriteRowNames',true);

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
        T10 = table(Multipliers_2x2, Adders_2x2, Multipliers_3x3, Pre_sum_3x3, Sum_3x3, Multipliers_4x4, Pre_sum_4x4, Sum4x4, Mult5x5, Pre_sum1_5x5, Pre_sum2_5x5, Det_5x5, 'RowNames', num_determinante);
        writetable(T10,['src/width_txt/Максимальные_значения_определителя_АЦП_№' num2str(i) '.xlsx'],'WriteRowNames',true)  


        for j = 1:length(num_determinante)
            width_mult_det2x2(j) = define_of_width_int(Multipliers_2x2(j), sim_options.type_2x2_det, sim_options.width_fractional);
            width_sum_det2x2(j) = define_of_width_int(Adders_2x2(j), sim_options.type_2x2_det, sim_options.width_fractional);

            width_mult_det3x3(j) = define_of_width_int(Multipliers_3x3(j), sim_options.type_3x3_det, sim_options.width_hilbert);
            width_presum_det3x3(j) = define_of_width_int(Pre_sum_3x3(j), sim_options.type_3x3_det, sim_options.width_hilbert);
            width_sum_det3x3(j) = define_of_width_int(Sum_3x3(j), sim_options.type_3x3_det, sim_options.width_hilbert);

            width_mult_det4x4(j) = define_of_width_int(Multipliers_4x4(j), sim_options.type_4x4_det, sim_options.width_hilbert);
            width_presum_det4x4(j) = define_of_width_int(Pre_sum_4x4(j), sim_options.type_4x4_det, sim_options.width_hilbert);
            width_sum_det4x4(j) = define_of_width_int(Sum4x4(j), sim_options.type_4x4_det, sim_options.width_hilbert);

            width_mult_det5x5(j) = define_of_width_int(Mult5x5(j), sim_options.type_5x5_det, sim_options.width_hilbert);
            width_presum1_det5x5(j) = define_of_width_int(Pre_sum1_5x5(j), sim_options.type_5x5_det, sim_options.width_hilbert);
            width_presum2_det5x5(j) = define_of_width_int(Pre_sum2_5x5(j), sim_options.type_5x5_det, sim_options.width_hilbert);
            width_sum_det5x5(j) = define_of_width_int(Det_5x5(j), sim_options.type_5x5_det, sim_options.width_hilbert);
        end

        % формируем таблицу разрядностей определителя 
        T11 = table(width_mult_det2x2, width_sum_det2x2, width_mult_det3x3, width_presum_det3x3, width_sum_det3x3, width_mult_det4x4, width_presum_det4x4, ...
            width_sum_det4x4, width_mult_det5x5, width_presum1_det5x5, width_presum2_det5x5, width_sum_det5x5, 'RowNames', num_determinante);
        writetable(T11,['src/width_txt/Разрядность_максимальных_значений_определителя_АЦП_№' num2str(i) '.xlsx'],'WriteRowNames',true)  

        %% Запись данных для делителя
        Dividend = Det_5x5(1,1);
        Divisor = Det_x3_int_max_in_cycle(:,i);
        Quotient = Divide_max_in_cycle(:,i);
        % формируем таблицу делителя
        T12 = table(Dividend, Divisor, Quotient, 'RowNames', num_divide);
        writetable(T12,['src/width_txt/Максимальные_значения_делителя_АЦП_№' num2str(i) '.xlsx'],'WriteRowNames',true)  
        
        %% Запись данных для адаптивного фильтра

        Multipliers_adaptive = Adaptive_filter_mult_array_max_in_cycle(:,i);
        Sum1_adaptive = Adaptive_filter_sum_array_max_in_cycle(:,1,i);
        Sum2_adaptive = Adaptive_filter_sum_array_max_in_cycle(:,2,i);
        Sum3_adaptive = Adaptive_filter_sum_array_max_in_cycle(:,3,i);

        % формируем таблицу адаптивного фильтра
        T13 = table(Multipliers_adaptive, Sum1_adaptive, Sum2_adaptive, Sum3_adaptive, 'RowNames', num_adaptive);
        writetable(T13,['src/width_txt/Максимальные_значения_адаптивного_фильтра_АЦП_№' num2str(i) '.xlsx'],'WriteRowNames',true); 

        Multipliers_adaptive_total_width = Adaptive_filter_mult_total_width_in_cycle(:,i);
        Sum1_adaptive_total_width = Adaptive_filter_sum_total_width_in_cycle(:,1,i);
        Sum2_adaptive_total_width = Adaptive_filter_sum_total_width_in_cycle(:,2,i);
        Sum3_adaptive_total_width = Adaptive_filter_sum_total_width_in_cycle(:,3,i);

        % формируем таблицу адаптивного фильтра
        T14 = table(Multipliers_adaptive_total_width, Sum1_adaptive_total_width, Sum2_adaptive_total_width, Sum3_adaptive_total_width, 'RowNames', num_adaptive);
        writetable(T14,['src/width_txt/Разрядность_элементов_адаптивного_фильтра_АЦП_№' num2str(i) '.xlsx'],'WriteRowNames',true);


        for j = 1:length(num_adaptive)
            width_mult_adaptive(j) = define_of_width_int(Multipliers_adaptive(j), sim_options.type_mult_in_adaptive_filter, sim_options.width_hilbert);
            width_sum1_adaptive(j) = define_of_width_int(Sum1_adaptive(j), sim_options.type_add_in_adaptive_filter, sim_options.width_hilbert);
            width_sum2_adaptive(j) = define_of_width_int(Sum2_adaptive(j), sim_options.type_add_in_adaptive_filter, sim_options.width_hilbert);
            width_sum3_adaptive(j) = define_of_width_int(Sum3_adaptive(j), sim_options.type_add_in_adaptive_filter, sim_options.width_hilbert);
        end

        % формируем таблицу разрядностей адаптивного фильтра
        T15 = table(width_mult_adaptive, width_sum1_adaptive, width_sum2_adaptive, width_sum3_adaptive, 'RowNames', num_adaptive);
        writetable(T15,['src/width_txt/Разрядность_максимальных_значений_адаптивного_фильтра_АЦП_№' num2str(i) '.xlsx'],'WriteRowNames',true)  

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