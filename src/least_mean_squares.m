function [y_array, y_array_double, y_array_int, DetM_2x2_array, DetM_2x2_array_int, ...
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
    ] = least_mean_square(adc_input, yri_cut, yri_cut_int, adaptive_mult_file, adaptive_sum_file, sim_options)

bb = 0;
vv = 0;
kk = 0;
y_outd = 0;

nn = 0;
x3 = cast(zeros(sim_options.Size_matrix,sim_options.Size_matrix), "double");                     
x3_int = cast(zeros(sim_options.Size_matrix,sim_options.Size_matrix), sim_options.type_fir_out); 

DetM_2x2_multiplier_total_abs_max = cast(zeros(sim_options.num_det2x2*2,1), sim_options.type_2x2_det);
Det2x2_sum_abs_max = cast(zeros(sim_options.num_det2x2,1), sim_options.type_2x2_det);

Mult_DetM_3x3_array_max = cast(zeros(sim_options.num_det2x2*3,1), sim_options.type_3x3_det);
Mult_DetM_3x3_array_mult_total_width_max = cast(zeros(sim_options.num_det2x2*3,1), sim_options.type_3x3_det);
DetM_3x3_int_pre_sum_array_max = cast(zeros(sim_options.num_det2x2,1), sim_options.type_3x3_det);
DetM_3x3_int_pre_sum_width_total_max = cast(zeros(sim_options.num_det2x2,1), sim_options.type_3x3_det);
DetM_3x3_int_sum_array_max = cast(zeros(sim_options.num_det2x2,1), sim_options.type_3x3_det);
DetM_3x3_int_sum_width_total_max = cast(zeros(sim_options.num_det2x2,1), sim_options.type_3x3_det);
% умножители определителя 4х4
DetM_4x4_int_mult_array_max = cast(zeros(sim_options.num_det2x2*2,1), sim_options.type_4x4_det);
% разрядность умножителей определителя 4х4
DetM_4x4_int_mult_width_total_max = cast(zeros(sim_options.num_det2x2*2,1), sim_options.type_4x4_det);
% пресумматоры определителя 4х4
DetM_4x4_int_pre_sum_array_max = cast(zeros(sim_options.num_det2x2,1), sim_options.type_4x4_det);
% разрядность пресумматоров определителя 4х4
DetM_4x4_int_pre_sum_width_total_max = cast(zeros(sim_options.num_det2x2,1), sim_options.type_4x4_det);
% сумматоры определителя 4х4
DetM_4x4_int_sum_array_max = cast(zeros(sim_options.Size_matrix,1), sim_options.type_4x4_det);
% разрядность сумматоров определителя 4х4
DetM_4x4_int_sum_array_width_total_max = cast(zeros(sim_options.Size_matrix,1), sim_options.type_4x4_det);
% умножители определителя 5х5
DetM_5x5_int_mult_array_max = cast(zeros(sim_options.Size_matrix,1), sim_options.type_5x5_det);
% разрядность умножителей 5x5
DetM_5x5_int_mult_array_width_total_max = cast(zeros(sim_options.Size_matrix,1), sim_options.type_5x5_det);
% пресумматоры1 определителя 5х5
DetM_5x5_int_pre_sum1_array_max = cast(zeros(2,1), sim_options.type_5x5_det);
% разрядность пресумматоров1 определителя 5х5
DetM_5x5_int_pre_sum1_array_width_total_max = cast(zeros(2,1), sim_options.type_5x5_det);
% пресумматор2 определителя 5х5
DetM_5x5_int_sum3_abs_max = cast(0, sim_options.type_5x5_det);
% разрядность пресумматора2 определителя 5х5
DetM_5x5_int_sum3_width_total_max = cast(0, sim_options.type_5x5_det);
% сумматор определителя 5х5
DetM_5x5_int_abs_max = cast(0, sim_options.type_5x5_det);
% разрядность сумматора определителя 5х5
DetM_5x5_int_width_total_max = cast(0, sim_options.type_5x5_det);
% начальный определитель
Det_x3_int_max = cast(0, sim_options.type_5x5_det);
% выход делителя
Divide_max = cast(0, sim_options.type_divide_out);

www1 = zeros(sim_options.Size_matrix,1);
www1_double = zeros(sim_options.Size_matrix,1);
www1_double_abs = zeros(sim_options.Size_matrix,1);
www1_int = cast(zeros(sim_options.Size_matrix,1), sim_options.type_divide_out);
www1_int_abs = cast(zeros(sim_options.Size_matrix,1), sim_options.type_divide_out);

Adaptive_filter_mult_array_max = cast(zeros(sim_options.Size_matrix,1),sim_options.type_mult_in_adaptive_filter);
Adaptive_filter_mult_total_width = cast(zeros(sim_options.Size_matrix,1),sim_options.type_mult_in_adaptive_filter);
Adaptive_filter_sum_array_max = cast(zeros(sim_options.Size_matrix,1),sim_options.type_add_in_adaptive_filter);
Adaptive_filter_sum_total_width = cast(zeros(sim_options.Size_matrix,1),sim_options.type_add_in_adaptive_filter);

for j = 1:sim_options.Size_matrix:length(yri_cut(:,1))-sim_options.Size_matrix+1 

    x3(1,:) = adc_input(j:sim_options.Size_matrix+j-1);
    x3(2,:) = adc_input(j+1:sim_options.Size_matrix+j);
    x3(3,:) = adc_input(j+2:sim_options.Size_matrix+j+1);
    x3(4,:) = adc_input(j+3:sim_options.Size_matrix+j+2);
    x3(5,:) = adc_input(j+4:sim_options.Size_matrix+j+3);

    x3_int(1,:) = adc_input(j:sim_options.Size_matrix+j-1);
    x3_int(2,:) = adc_input(j+1:sim_options.Size_matrix+j);
    x3_int(3,:) = adc_input(j+2:sim_options.Size_matrix+j+1);
    x3_int(4,:) = adc_input(j+3:sim_options.Size_matrix+j+2);
    x3_int(5,:) = adc_input(j+4:sim_options.Size_matrix+j+3);

    %% determinant    
    det_matlab(j) = det(x3);

	[det_x3(j), det_x3_int(j), Det_5x5_LU_matlab_array(j), DetM_2x2, DetM_2x2_int, Det_2x2_LU_matlab, DetM_3x3_array, Det_3x3_LU_matlab, DetM_4x4_array, Det_4x4_LU_matlab, ...
		... % умножители определителя 2x2 
		DetM_2x2_multiplier_total_abs, ...
		... % сумматоры определителя 2x2 
		Det2x2_sum_abs, ...
		... % умножители определителя 3х3
		Mult_DetM_3x3_array, ...
		... % разрядность умножителей определителя 3x3
		Mult_DetM_3x3_array_mult_total_width, ...
		... % пресумматоры определителя 3х3
		DetM_3x3_int_pre_sum_array, ...
		... % разрядность пресумматоров определителя 3x3
		DetM_3x3_int_pre_sum_width_total, ...
		... % сумматоры определителя 3х3
		DetM_3x3_int_sum_array, ...
		... % разрядность сумматоров определителя 3x3
		DetM_3x3_int_sum_width_total, ...
		... % умножители определителя 4х4
		DetM_4x4_int_mult_array, ...
		... % разрядность умножителей определителя 4х4
		DetM_4x4_int_mult_width_total, ...
		... % пресумматоры определителя 4х4
		DetM_4x4_int_pre_sum_array, ...
		... % разрядность пресумматоров определителя 4х4
		DetM_4x4_int_pre_sum_width_total, ...
		... % сумматоры определителя 4х4
		DetM_4x4_int_sum_array, ...
		... % разрядность сумматоров определителя 4х4
		DetM_4x4_int_sum_array_width_total, ...
		... % умножители определителя 5х5
		DetM_5x5_int_mult_array, ...
		... % разрядность умножителей 5x5
		DetM_5x5_int_mult_array_width_total, ...
		... % пресумматоры1 определителя 5х5
		DetM_5x5_int_pre_sum1_array, ...
		... % разрядность пресумматоров1 определителя 5х5
		DetM_5x5_int_pre_sum1_array_width_total, ...
		... % пресумматор2 определителя 5х5
		DetM_5x5_int_sum3_abs, ...
		... % разрядность пресумматора2 определителя 5х5
		DetM_5x5_int_sum3_width_total, ...
		... % сумматор определителя 5х5
		DetM_5x5_int_abs, ...
		... % разрядность сумматора определителя 5х5
		DetM_5x5_int_width_total ...
	] = determinate(x3, x3_int, sim_options); % int

    if (det_matlab(j) == 0)
        det_matlab(j) = 1;
        disp('Determinant LU x3 equal 0');
    end
    if (det_x3_int(j) == 0)
        det_x3_int(j) = 1;
        disp('Determinant int x3 equal 0');
    end
    if (det_x3(j) == 0)
        det_x3(j) = 1;
        disp('Determinant double x3 equal 0');
    end

	DetM_2x2_array(bb+1:bb+10) = abs(DetM_2x2);
	DetM_2x2_array_int(bb+1:bb+10) = DetM_2x2_int;
	Det_2x2_LU_matlab_array(bb+1:bb+10) = Det_2x2_LU_matlab;
	
	DetM_3x3_array_int(bb+1:bb+10) = DetM_3x3_int_sum_array;
	DetM_3x3_array_dd(bb+1:bb+10) = DetM_3x3_array;
	Det_3x3_LU_matlab_array(bb+1:bb+10) = Det_3x3_LU_matlab;

	DetM_4x4_array_int(bb+1:bb+5) = DetM_4x4_int_sum_array;
	DetM_4x4_array_dd(bb+1:bb+5) = DetM_4x4_array;
	Det_4x4_LU_matlab_array(bb+1:bb+5) = Det_4x4_LU_matlab;

	bb = bb + 10;
	vv = vv + 5;

	%% определяем макс. значения в определителе 2x2
	for n = 1:sim_options.num_det2x2*2
		% определяем максимальное значение на каждом из 20 умножителей
		if DetM_2x2_multiplier_total_abs_max(n) < DetM_2x2_multiplier_total_abs(n) 
			DetM_2x2_multiplier_total_abs_max(n) = DetM_2x2_multiplier_total_abs(n);
		end
	end 
	for n = 1:sim_options.num_det2x2
		% определяем максимальное значение на каждом из 10 сумматоров
		if Det2x2_sum_abs_max(n) < Det2x2_sum_abs(n) 
			Det2x2_sum_abs_max(n) = Det2x2_sum_abs(n);
		end
	end 
	%% 																		3x3
	for n = 1:sim_options.num_det2x2*3
		% определяем максимальное значение на каждом из 30 умножителей
		if Mult_DetM_3x3_array_max(n) < Mult_DetM_3x3_array(n) 
			Mult_DetM_3x3_array_max(n) = Mult_DetM_3x3_array(n);
		end		
		% определяем макс. разрядность на каждом из 30 умножителей
		if Mult_DetM_3x3_array_mult_total_width_max(n) < Mult_DetM_3x3_array_mult_total_width(n) 
			Mult_DetM_3x3_array_mult_total_width_max(n) = Mult_DetM_3x3_array_mult_total_width(n);
		end
	end
	for n = 1:sim_options.num_det2x2
		% определяем максимальное значение на каждом из 10 пресумматорах
		if DetM_3x3_int_pre_sum_array_max(n) < DetM_3x3_int_pre_sum_array(n) 
			DetM_3x3_int_pre_sum_array_max(n) = DetM_3x3_int_pre_sum_array(n);
		end		
		% определяем макс. разрядность на каждом из 10 пресумматорах
		if DetM_3x3_int_pre_sum_width_total_max(n) < DetM_3x3_int_pre_sum_width_total(n) 
			DetM_3x3_int_pre_sum_width_total_max(n) = DetM_3x3_int_pre_sum_width_total(n);
		end
		% определяем максимальное значение на каждом из 10 сумматорах
		if DetM_3x3_int_sum_array_max(n) < DetM_3x3_int_sum_array(n) 
			DetM_3x3_int_sum_array_max(n) = DetM_3x3_int_sum_array(n);
		end		
		% определяем макс. разрядность на каждом из 10 сумматорах
		if DetM_3x3_int_sum_width_total_max(n) < DetM_3x3_int_sum_width_total(n) 
			DetM_3x3_int_sum_width_total_max(n) = DetM_3x3_int_sum_width_total(n);
		end
	end
	%% 																		4x4 
	for n = 1:sim_options.num_det2x2*2
		% определяем максимальное значение на каждом из 20 умножителей
		if  DetM_4x4_int_mult_array_max(n) < DetM_4x4_int_mult_array(n) 
			DetM_4x4_int_mult_array_max(n) = DetM_4x4_int_mult_array(n);
		end		
		% определяем макс. разрядность на каждом из 20 умножителей
		if  DetM_4x4_int_mult_width_total_max(n) < DetM_4x4_int_mult_width_total(n) 
			DetM_4x4_int_mult_width_total_max(n) = DetM_4x4_int_mult_width_total(n);
		end
	end
	for n = 1:sim_options.num_det2x2
		% определяем максимальное значение на каждом из 10 пресумматоров
		if  DetM_4x4_int_pre_sum_array_max(n) < DetM_4x4_int_pre_sum_array(n) 
			DetM_4x4_int_pre_sum_array_max(n) = DetM_4x4_int_pre_sum_array(n);
		end		
		% определяем макс. разрядность на каждом из 10 пресумматоров
		if  DetM_4x4_int_pre_sum_width_total_max(n) < DetM_4x4_int_pre_sum_width_total(n) 
			DetM_4x4_int_pre_sum_width_total_max(n) = DetM_4x4_int_pre_sum_width_total(n);
		end
	end
	for n = 1:5
		% определяем максимальное значение на каждом из 5 сумматоров
		if  DetM_4x4_int_sum_array_max(n) < DetM_4x4_int_sum_array(n) 
			DetM_4x4_int_sum_array_max(n) = DetM_4x4_int_sum_array(n);
		end		
		% определяем макс. разрядность на каждом из 5 сумматоров
		if  DetM_4x4_int_sum_array_width_total(n) < DetM_4x4_int_sum_array_width_total(n) 
			DetM_4x4_int_sum_array_width_total(n) = DetM_4x4_int_sum_array_width_total(n);
		end
	end
	%% 																		5x5
	for n = 1:5
		% определяем максимальное значение на каждом из 5 сумматоров
		if  DetM_5x5_int_mult_array_max(n) < DetM_5x5_int_mult_array(n) 
			DetM_5x5_int_mult_array_max(n) = DetM_5x5_int_mult_array(n);
		end		
		% определяем макс. разрядность на каждом из 5 сумматоров
		if  DetM_5x5_int_mult_array_width_total_max(n) < DetM_5x5_int_mult_array_width_total(n) 
			DetM_5x5_int_mult_array_width_total_max(n) = DetM_5x5_int_mult_array_width_total(n);
		end
	end
	for n = 1:2
		% определяем максимальное значение на каждом из 5 сумматоров
		if  DetM_5x5_int_pre_sum1_array_max(n) < DetM_5x5_int_pre_sum1_array(n) 
			DetM_5x5_int_pre_sum1_array_max(n) = DetM_5x5_int_pre_sum1_array(n);
		end		
		% определяем макс. разрядность на каждом из 5 сумматоров
		if  DetM_5x5_int_pre_sum1_array_width_total_max(n) < DetM_5x5_int_pre_sum1_array_width_total(n) 
			DetM_5x5_int_pre_sum1_array_width_total_max(n) = DetM_5x5_int_pre_sum1_array_width_total(n);
		end
	end
	% пресумматор2 определителя 5х5
	if  DetM_5x5_int_sum3_abs_max < DetM_5x5_int_sum3_abs 
		DetM_5x5_int_sum3_abs_max = DetM_5x5_int_sum3_abs;
	end
	% разрядность пресумматора2 определителя 5х5
	if  DetM_5x5_int_sum3_width_total_max < DetM_5x5_int_sum3_width_total 
		DetM_5x5_int_sum3_width_total_max = DetM_5x5_int_sum3_width_total;
	end
	% сумматор определителя 5х5
	if  DetM_5x5_int_abs_max < DetM_5x5_int_abs 
		DetM_5x5_int_abs_max = DetM_5x5_int_abs;
	end
	% разрядность сумматора определителя 5х5
	if  DetM_5x5_int_width_total_max < DetM_5x5_int_width_total 
		DetM_5x5_int_width_total_max = DetM_5x5_int_width_total;
    end
    % начальный определитель
    if Det_x3_int_max < DetM_5x5_int_abs
        Det_x3_int_max = DetM_5x5_int_abs;
    end
	%% 

	for i = 1:sim_options.Size_matrix
		kk = kk + 1;

		x3_shift = x3;
		x3_shift(1:sim_options.Size_matrix,i) = yri_cut(j:j+sim_options.Size_matrix-1);

		x3_shift_int = x3_int;
		x3_shift_int(1:sim_options.Size_matrix,i) = yri_cut_int(j:j+sim_options.Size_matrix-1); 

		w1 = lsqminnorm(x3, yri_cut(j:sim_options.Size_matrix+j-1));
		%% determinant
		det_x3_shift(kk) = det(x3_shift);
		[det_out_shift(kk), det_out_shift_int(kk), Det_5x5_LU_matlab_array(j), DetM_2x2, DetM_2x2_int, Det_2x2_LU_matlab, DetM_3x3_array, Det_3x3_LU_matlab, DetM_4x4_array, Det_4x4_LU_matlab, ...
			... % умножители определителя 2x2 
			DetM_2x2_multiplier_total_abs, ...
			... % сумматоры определителя 2x2 
			Det2x2_sum_abs, ...
			... % умножители определителя 3х3
			Mult_DetM_3x3_array, ...
			... % разрядность умножителей определителя 3x3
			Mult_DetM_3x3_array_mult_total_width, ...
			... % пресумматоры определителя 3х3
			DetM_3x3_int_pre_sum_array, ...
			... % разрядность пресумматоров определителя 3x3
			DetM_3x3_int_pre_sum_width_total, ...
			... % сумматоры определителя 3х3
			DetM_3x3_int_sum_array, ...
			... % разрядность сумматоров определителя 3x3
			DetM_3x3_int_sum_width_total, ...
			... % умножители определителя 4х4
			DetM_4x4_int_mult_array, ...
			... % разрядность умножителей определителя 4х4
			DetM_4x4_int_mult_width_total, ...
			... % пресумматоры определителя 4х4
			DetM_4x4_int_pre_sum_array, ...
			... % разрядность пресумматоров определителя 4х4
			DetM_4x4_int_pre_sum_width_total, ...
			... % сумматоры определителя 4х4
			DetM_4x4_int_sum_array, ...
			... % разрядность сумматоров определителя 4х4
			DetM_4x4_int_sum_array_width_total, ...
			... % умножители определителя 5х5
			DetM_5x5_int_mult_array, ...
			... % разрядность умножителей 5x5
			DetM_5x5_int_mult_array_width_total, ...
			... % пресумматоры1 определителя 5х5
			DetM_5x5_int_pre_sum1_array, ...
			... % разрядность пресумматоров1 определителя 5х5
			DetM_5x5_int_pre_sum1_array_width_total, ...
			... % пресумматор2 определителя 5х5
			DetM_5x5_int_sum3_abs, ...
			... % разрядность пресумматора2 определителя 5х5
			DetM_5x5_int_sum3_width_total, ...
			... % сумматор определителя 5х5
			DetM_5x5_int_abs, ...
			... % разрядность сумматора определителя 5х5
			DetM_5x5_int_width_total ...
		] = determinate(x3_shift, x3_shift_int, sim_options);

        if (det_x3_shift(kk) == 0)
            det_x3_shift(kk) = 1;
            disp('Determinant LU x3_shift equal 0');
        end
        if (det_out_shift_int(kk) == 0)
            det_out_shift_int(kk) = 1;
            disp('Determinant int x3_shift equal 0');
        end
        if (det_out_shift(kk) == 0)
            det_out_shift(kk) = 1;
            disp('Determinant double x3_shift equal 0');
        end

        %%
		if i == 1
		   DetM_2x2_int_double = [cast(DetM_2x2_int(1:4), sim_options.type_2x2_det);  cast(DetM_2x2_int(5:10), sim_options.type_2x2_det)];
		   % DetM_2x2_double = [(DetM_2x2(1:4))*2^-remainder;  (DetM_2x2(5:10))];

		   DetM_3x3_int_double = [cast(DetM_3x3_int_sum_array(1:4),sim_options.type_3x3_det); cast(DetM_3x3_int_sum_array(5:10),sim_options.type_3x3_det)];

		   DetM_4x4_int_double = [cast(DetM_4x4_int_sum_array(1),sim_options.type_4x4_det); cast(DetM_4x4_int_sum_array(2:5),sim_options.type_4x4_det)];
		elseif i == 2
		   DetM_2x2_int_double = [cast(DetM_2x2_int(1),sim_options.type_2x2_det); cast(DetM_2x2_int(2:4),sim_options.type_2x2_det); cast(DetM_2x2_int(5:7),sim_options.type_2x2_det); cast(DetM_2x2_int(8:10),sim_options.type_2x2_det)];
		   DetM_2x2_double = [(DetM_2x2(1)); (DetM_2x2(2:4)); (DetM_2x2(5:7)); (DetM_2x2(8:10))];

		   DetM_3x3_int_double = [cast(DetM_3x3_int_sum_array(1),sim_options.type_3x3_det); cast(DetM_3x3_int_sum_array(2:4),sim_options.type_3x3_det); cast(DetM_3x3_int_sum_array(5:7),sim_options.type_3x3_det); cast(DetM_3x3_int_sum_array(8:10),sim_options.type_3x3_det)];

		   DetM_4x4_int_double = [cast(DetM_4x4_int_sum_array(1),sim_options.type_4x4_det); cast(DetM_4x4_int_sum_array(2),sim_options.type_4x4_det); cast(DetM_4x4_int_sum_array(3:5),sim_options.type_4x4_det)];
		elseif i == 3
		   DetM_2x2_int_double = [cast(DetM_2x2_int(1),sim_options.type_2x2_det); cast(DetM_2x2_int(2),sim_options.type_2x2_det); cast(DetM_2x2_int(3:4),sim_options.type_2x2_det); cast(DetM_2x2_int(5),sim_options.type_2x2_det); cast(DetM_2x2_int(6:7),sim_options.type_2x2_det); ...
		   cast(DetM_2x2_int(8:9),sim_options.type_2x2_det); cast(DetM_2x2_int(10),sim_options.type_2x2_det)];

		   % DetM_2x2_double = [(DetM_2x2(1)); (DetM_2x2(2))*2^-remainder; (DetM_2x2(3:4)); (DetM_2x2(5))*2^-remainder; (DetM_2x2(6:7)); ...
		   % cast(DetM_2x2(8:9),sim_options.type_2x2_det)*2^-remainder; (DetM_2x2(10))]; 

		   DetM_3x3_int_double = [cast(DetM_3x3_int_sum_array(1),sim_options.type_3x3_det); cast(DetM_3x3_int_sum_array(2),sim_options.type_3x3_det); cast(DetM_3x3_int_sum_array(3:4),sim_options.type_3x3_det); cast(DetM_3x3_int_sum_array(5),sim_options.type_3x3_det); 
		       cast(DetM_3x3_int_sum_array(6:7),sim_options.type_3x3_det); cast(DetM_3x3_int_sum_array(8:9),sim_options.type_3x3_det); cast(DetM_3x3_int_sum_array(10),sim_options.type_3x3_det)];

		   DetM_4x4_int_double = [cast(DetM_4x4_int_sum_array(1:2),sim_options.type_4x4_det); cast(DetM_4x4_int_sum_array(3),sim_options.type_4x4_det); cast(DetM_4x4_int_sum_array(4:5),sim_options.type_4x4_det)];
		elseif i == 4
		   DetM_2x2_int_double = [cast(DetM_2x2_int(1:2),sim_options.type_2x2_det); cast(DetM_2x2_int(3),sim_options.type_2x2_det); cast(DetM_2x2_int(4:5),sim_options.type_2x2_det); cast(DetM_2x2_int(6),sim_options.type_2x2_det); cast(DetM_2x2_int(7),sim_options.type_2x2_det); ...
		   cast(DetM_2x2_int(8),sim_options.type_2x2_det); cast(DetM_2x2_int(9),sim_options.type_2x2_det); cast(DetM_2x2_int(10),sim_options.type_2x2_det)];

		   % DetM_2x2_double = [(DetM_2x2(1:2)); (DetM_2x2(3))*2^-remainder; (DetM_2x2(4:5)); (DetM_2x2(6))*2^-remainder; (DetM_2x2(7)); ...
		   % (DetM_2x2(8))*2^-remainder; (DetM_2x2(9)); (DetM_2x2(10))*2^-remainder];

		   DetM_3x3_int_double = [cast(DetM_3x3_int_sum_array(1:2),sim_options.type_3x3_det); cast(DetM_3x3_int_sum_array(3),sim_options.type_3x3_det); cast(DetM_3x3_int_sum_array(4:5),sim_options.type_3x3_det); cast(DetM_3x3_int_sum_array(6),sim_options.type_3x3_det); cast(DetM_3x3_int_sum_array(7),sim_options.type_3x3_det); ...
		   cast(DetM_3x3_int_sum_array(8),sim_options.type_3x3_det); cast(DetM_3x3_int_sum_array(9),sim_options.type_3x3_det); cast(DetM_3x3_int_sum_array(10),sim_options.type_3x3_det)];

		   DetM_4x4_int_double = [cast(DetM_4x4_int_sum_array(1:3),sim_options.type_4x4_det); cast(DetM_4x4_int_sum_array(4),sim_options.type_4x4_det); cast(DetM_4x4_int_sum_array(5),sim_options.type_4x4_det)];
		elseif i == 5
		   DetM_2x2_int_double = [cast(DetM_2x2_int(1:3),sim_options.type_2x2_det); cast(DetM_2x2_int(4),sim_options.type_2x2_det); cast(DetM_2x2_int(5:6),sim_options.type_2x2_det); cast(DetM_2x2_int(7),sim_options.type_2x2_det); cast(DetM_2x2_int(8),sim_options.type_2x2_det); ...
		   cast(DetM_2x2_int(9:10),sim_options.type_2x2_det)];

		   % DetM_2x2_double = [(DetM_2x2(1:3)); (DetM_2x2(4))*2^-remainder; (DetM_2x2(5:6)); (DetM_2x2(7))*2^-remainder; (DetM_2x2(8)); ...
		   % (DetM_2x2(9:10))*2^-remainder];

		   DetM_3x3_int_double = [cast(DetM_3x3_int_sum_array(1:3),sim_options.type_3x3_det); cast(DetM_3x3_int_sum_array(4),sim_options.type_3x3_det); cast(DetM_3x3_int_sum_array(5:6),sim_options.type_3x3_det); cast(DetM_3x3_int_sum_array(7),sim_options.type_3x3_det); cast(DetM_3x3_int_sum_array(8),sim_options.type_3x3_det); ...
		   cast(DetM_3x3_int_sum_array(9:10),sim_options.type_3x3_det)];

		   DetM_4x4_int_double = [cast(DetM_4x4_int_sum_array(1:4),sim_options.type_4x4_det); cast(DetM_4x4_int_sum_array(5),sim_options.type_4x4_det)];
		end

		DetM_2x2_array(bb+1:bb+10) = abs(DetM_2x2);
		DetM_2x2_array_int(bb+1:bb+10) = DetM_2x2_int_double;
		Det_2x2_LU_matlab_array(bb+1:bb+10) = Det_2x2_LU_matlab;
	
		DetM_3x3_array_int(bb+1:bb+10) = DetM_3x3_int_double;
		DetM_3x3_array_dd(bb+1:bb+10) = DetM_3x3_array;
		Det_3x3_LU_matlab_array(bb+1:bb+10) = Det_3x3_LU_matlab;
	
		DetM_4x4_array_int(bb+1:bb+5) = DetM_4x4_int_double;
		DetM_4x4_array_dd(bb+1:bb+5) = DetM_4x4_array;
		Det_4x4_LU_matlab_array(bb+1:bb+5) = Det_4x4_LU_matlab;

		bb = bb + 10;
		vv = vv + 5;

		%% определяем макс. значения в определителе 2x2
		for n = 1:sim_options.num_det2x2*2
			% определяем максимальное значение на каждом из 20 умножителей
			if DetM_2x2_multiplier_total_abs_max(n) < DetM_2x2_multiplier_total_abs(n) 
				DetM_2x2_multiplier_total_abs_max(n) = DetM_2x2_multiplier_total_abs(n);
			end
		end 
		for n = 1:sim_options.num_det2x2
			% определяем максимальное значение на каждом из 10 сумматоров
			if Det2x2_sum_abs_max(n) < Det2x2_sum_abs(n) 
				Det2x2_sum_abs_max(n) = Det2x2_sum_abs(n);
			end
		end 
		%% 																		3x3
		for n = 1:sim_options.num_det2x2*3
			% определяем максимальное значение на каждом из 30 умножителей
			if Mult_DetM_3x3_array_max(n) < Mult_DetM_3x3_array(n) 
				Mult_DetM_3x3_array_max(n) = Mult_DetM_3x3_array(n);
			end		
			% определяем макс. разрядность на каждом из 30 умножителей
			if Mult_DetM_3x3_array_mult_total_width_max(n) < Mult_DetM_3x3_array_mult_total_width(n) 
				Mult_DetM_3x3_array_mult_total_width_max(n) = Mult_DetM_3x3_array_mult_total_width(n);
            end
		end
		for n = 1:sim_options.num_det2x2
			% определяем максимальное значение на каждом из 10 пресумматорах
			if DetM_3x3_int_pre_sum_array_max(n) < DetM_3x3_int_pre_sum_array(n) 
				DetM_3x3_int_pre_sum_array_max(n) = DetM_3x3_int_pre_sum_array(n);
			end		
			% определяем макс. разрядность на каждом из 10 пресумматорах
			if DetM_3x3_int_pre_sum_width_total_max(n) < DetM_3x3_int_pre_sum_width_total(n) 
				DetM_3x3_int_pre_sum_width_total_max(n) = DetM_3x3_int_pre_sum_width_total(n);
			end
			% определяем максимальное значение на каждом из 10 сумматорах
			if DetM_3x3_int_sum_array_max(n) < DetM_3x3_int_sum_array(n) 
				DetM_3x3_int_sum_array_max(n) = DetM_3x3_int_sum_array(n);
			end		
			% определяем макс. разрядность на каждом из 10 сумматорах
			if DetM_3x3_int_sum_width_total_max(n) < DetM_3x3_int_sum_width_total(n) 
				DetM_3x3_int_sum_width_total_max(n) = DetM_3x3_int_sum_width_total(n);
			end
		end
		%% 																		4x4 
		for n = 1:sim_options.num_det2x2*2
			% определяем максимальное значение на каждом из 20 умножителей
			if  DetM_4x4_int_mult_array_max(n) < DetM_4x4_int_mult_array(n) 
				DetM_4x4_int_mult_array_max(n) = DetM_4x4_int_mult_array(n);
			end		
			% определяем макс. разрядность на каждом из 20 умножителей
			if  DetM_4x4_int_mult_width_total_max(n) < DetM_4x4_int_mult_width_total(n) 
				DetM_4x4_int_mult_width_total_max(n) = DetM_4x4_int_mult_width_total(n);
			end
		end
		for n = 1:sim_options.num_det2x2
			% определяем максимальное значение на каждом из 10 пресумматоров
			if  DetM_4x4_int_pre_sum_array_max(n) < DetM_4x4_int_pre_sum_array(n) 
				DetM_4x4_int_pre_sum_array_max(n) = DetM_4x4_int_pre_sum_array(n);
			end		
			% определяем макс. разрядность на каждом из 10 пресумматоров
			if  DetM_4x4_int_pre_sum_width_total_max(n) < DetM_4x4_int_pre_sum_width_total(n) 
				DetM_4x4_int_pre_sum_width_total_max(n) = DetM_4x4_int_pre_sum_width_total(n);
			end
		end
		for n = 1:5
			% определяем максимальное значение на каждом из 5 сумматоров
			if  DetM_4x4_int_sum_array_max(n) < DetM_4x4_int_sum_array(n) 
				DetM_4x4_int_sum_array_max(n) = DetM_4x4_int_sum_array(n);
			end		
			% определяем макс. разрядность на каждом из 5 сумматоров
			if  DetM_4x4_int_sum_array_width_total(n) < DetM_4x4_int_sum_array_width_total(n) 
				DetM_4x4_int_sum_array_width_total(n) = DetM_4x4_int_sum_array_width_total(n);
			end
		end
		%% 																		5x5
		for n = 1:5
			% определяем максимальное значение на каждом из 5 сумматоров
			if  DetM_5x5_int_mult_array_max(n) < DetM_5x5_int_mult_array(n) 
				DetM_5x5_int_mult_array_max(n) = DetM_5x5_int_mult_array(n);
			end		
			% определяем макс. разрядность на каждом из 5 сумматоров
			if  DetM_5x5_int_mult_array_width_total_max(n) < DetM_5x5_int_mult_array_width_total(n) 
				DetM_5x5_int_mult_array_width_total_max(n) = DetM_5x5_int_mult_array_width_total(n);
			end
		end
		for n = 1:2
			% определяем максимальное значение на каждом из 5 сумматоров
			if  DetM_5x5_int_pre_sum1_array_max(n) < DetM_5x5_int_pre_sum1_array(n) 
				DetM_5x5_int_pre_sum1_array_max(n) = DetM_5x5_int_pre_sum1_array(n);
			end		
			% определяем макс. разрядность на каждом из 5 сумматоров
			if  DetM_5x5_int_pre_sum1_array_width_total_max(n) < DetM_5x5_int_pre_sum1_array_width_total(n) 
				DetM_5x5_int_pre_sum1_array_width_total_max(n) = DetM_5x5_int_pre_sum1_array_width_total(n);
			end
		end
		% пресумматор2 определителя 5х5
		if  DetM_5x5_int_sum3_abs_max < DetM_5x5_int_sum3_abs 
			DetM_5x5_int_sum3_abs_max = DetM_5x5_int_sum3_abs;
		end
		% разрядность пресумматора2 определителя 5х5
		if  DetM_5x5_int_sum3_width_total_max < DetM_5x5_int_sum3_width_total 
			DetM_5x5_int_sum3_width_total_max = DetM_5x5_int_sum3_width_total;
		end
		% сумматор определителя 5х5
		if  DetM_5x5_int_abs_max < DetM_5x5_int_abs 
			DetM_5x5_int_abs_max = DetM_5x5_int_abs;
		end
		% разрядность сумматора определителя 5х5
		if  DetM_5x5_int_width_total_max < DetM_5x5_int_width_total 
			DetM_5x5_int_width_total_max = DetM_5x5_int_width_total;
        end

		%% divide determinant
        www1(i,:) = det_x3_shift(kk) ./ det_matlab(j); 
        % double
        [www1_double(i,:), overflow_divide_double(i,:), www1_double_abs(i,:), width_total_double(i,:)] = ...
            divide(det_out_shift(kk), det_x3(j), "double", "double", 64, 64, "double", sim_options.divide_factor);
        % int
		[www1_int(i,:), overflow_divide_int(i,:), www1_int_abs(i,:), width_total_int(i,:)] = ...
            divide(det_out_shift_int(kk), det_x3_int(j), sim_options.type_5x5_det, sim_options.type_5x5_det, 64, 64, sim_options.type_divide_out, sim_options.divide_factor); % int

        if (overflow_divide_int(i,:) == 1)
            disp('Переполнение делителя');
            disp({sim_options.SNR, sim_options.freq});
            disp({www1_int(i,:), det_out_shift_int(kk), int64(det_x3_int(j))});
        end

        % находим макс.значение выхода делителя
        if Divide_max < (www1_int_abs(i,:)) 
            Divide_max = www1_int_abs(i,:);
        end
    end

	%% adaptive filter
    y_outd = w1(1).*x3(:,1)+w1(2).*x3(:,2)+w1(3).*x3(:,3)+w1(4).*x3(:,4)+w1(5).*x3(:,5);

    x3_int_c = cast(x3_int, sim_options.type_mult_in_adaptive_filter);
    www1_int_c = cast(www1_int, sim_options.type_mult_in_adaptive_filter);

    [y_out, y_out_int, y_int_abs, y_int_total_width, sum_array_out, sum_int_width_total] ...
    = adaptive_filter(x3, www1_double, x3_int_c, www1_int_c, adaptive_mult_file, adaptive_sum_file, sim_options);

    % y_out_int = bitshift(y_out_int, -sim_options.divide_factor); % сдвигаем данные

    y_out_int_shift = y_out_int;
    y_out_double = round(y_out * 2^-sim_options.divide_factor);

    % округление значений после фильтра
    y_out_int_shift = round_int(y_out_int_shift, sim_options.divide_factor, sim_options.type_add_in_adaptive_filter);


	%% определяем макс. значения
	for n = 1:sim_options.Size_matrix
		% определяем максимальное значение на каждом умножителе
		if Adaptive_filter_mult_array_max(n) < y_int_abs(n) 
			Adaptive_filter_mult_array_max(n) = y_int_abs(n);
        end
		% определяем максимальную разрядность умножителей
		if Adaptive_filter_mult_total_width(n) < y_int_total_width(n) 
			Adaptive_filter_mult_total_width(n) = y_int_total_width(n);
        end
		% определяем максимальное значение сумматоров
		if Adaptive_filter_sum_array_max(n) < sum_array_out(n) 
			Adaptive_filter_sum_array_max(n) = sum_array_out(n);
        end
        % определяем максимальную разрядность сумматоров
		if Adaptive_filter_sum_total_width(n) < sum_int_width_total(n) 
			Adaptive_filter_sum_total_width(n) = sum_int_width_total(n);
        end
	end 

	y_array(nn+1:nn+sim_options.Size_matrix,:) = y_outd;
    y_array_double(nn+1:nn+sim_options.Size_matrix,:) = y_out_double;
	y_array_int(nn+1:nn+sim_options.Size_matrix,:) = y_out_int_shift;
    nn = nn + sim_options.Size_matrix;

end

    %% 2x2
    % relativeError_DetM_2x2_Myfunc_double_vs_Myfunc_int = DetM_2x2_array./double(DetM_2x2_array_int);
	% relativeError_DetM_2x2_Myfunc_double_vs_Matlab_LU = Det_2x2_LU_matlab_array./double(DetM_2x2_array);
    % 
    % figure(18)
	% subplot(2,1,1)
	% plot(relativeError_DetM_2x2_Myfunc_double_vs_Matlab_LU, '-o');
	% title('Относительная ошибка между определителями 2x2, найденных с помощью прямого нахождения Int vs Функции Матлаб')
    % ylabel('Величина ошибки') 
    % xlabel('Номер отсчета') 
	% subplot(2,1,2)
    % plot(relativeError_DetM_2x2_Myfunc_double_vs_Myfunc_int, '-o');
    % title('Относительная ошибка между определителями 2x2, найденных с помощью прямого нахождения Double vs Single')
    % ylabel('Величина ошибки') 
    % xlabel('Номер отсчета') 
    % 
	%% 3x3
    % relativeError_DetM_3x3_Myfunc_double_vs_Myfunc_int = DetM_3x3_array_dd./double(DetM_3x3_array_int);
	% relativeError_DetM_3x3_Myfunc_double_vs_Matlab_LU = Det_3x3_LU_matlab_array./double(DetM_3x3_array_int);
    % 
    % figure(19)
	% subplot(2,1,1)
    % plot(relativeError_DetM_3x3_Myfunc_double_vs_Matlab_LU, '-o');
    % title('Относительная ошибка между определителями 3x3, найденных с помощью прямого нахождения Int vs Функции Матлаб')
    % ylabel('Величина ошибки') 
    % xlabel('Номер отсчета') 
	% subplot(2,1,2)
    % plot(relativeError_DetM_3x3_Myfunc_double_vs_Myfunc_int, '-o');
    % title('Относительная ошибка между определителями 3x3, найденных с помощью прямого нахождения Double vs Single')
    % ylabel('Величина ошибки') 
    % xlabel('Номер отсчета') 
	
	%% 4x4
	% relativeError_DetM_4x4_Myfunc_double_vs_Myfunc_int = DetM_4x4_array_dd./double(DetM_4x4_array_int);
	% relativeError_DetM_4x4_Myfunc_double_vs_Matlab_LU = Det_4x4_LU_matlab_array./double(DetM_4x4_array_int);
    % 
    % figure(20)
	% subplot(2,1,1)
    % plot(relativeError_DetM_4x4_Myfunc_double_vs_Matlab_LU, '-o');
    % title('Относительная ошибка между определителями 4x4, найденных с помощью прямого нахождения Int(пока double) vs Функции Матлаб')
    % ylabel('Величина ошибки') 
    % xlabel('Номер отсчета') 
	% subplot(2,1,2)
    % plot(relativeError_DetM_4x4_Myfunc_double_vs_Myfunc_int, '-o');
    % title('Относительная ошибка между определителями 4x4, найденных с помощью прямого нахождения Int(пока double) vs Double')
    % ylabel('Величина ошибки') 
    % xlabel('Номер отсчета') 

    %% 5x5
	relativeError_DetM_5x5_LU_vs_Myfunc = det_x3_shift./double(det_out_shift_int);
    relativeError_DetM_5x5_Myfunc_double_vs_Myfunc_int = det_out_shift./double(det_out_shift_int);

	figure(21)
    subplot(2,1,1)
	plot(relativeError_DetM_5x5_LU_vs_Myfunc, '-o');
    title('Относительная ошибка между определителями, найденных с помощью LU-преобразования и прямого нахождения')
    ylabel('Величина ошибки') 
    xlabel('Номер отсчета') 
    subplot(2,1,2)
	plot(relativeError_DetM_5x5_Myfunc_double_vs_Myfunc_int, '-o');
	title('Относительная ошибка между определителями, найденных с помощью прямого нахождения Double vs Single')
    ylabel('Величина ошибки') 
    xlabel('Номер отсчета') 

    relativeError_y_out_double_vs_y_out_int = y_array_double ./ double(y_array_int);
    figure(29); plot(relativeError_y_out_double_vs_y_out_int); %*2^-(sim_options.divide_factor));

end