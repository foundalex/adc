function [y_array, y_array_int, DetM_2x2_array, DetM_2x2_array_int, ...
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
	DetM_5x5_int_width_total_max ...
    ] = least_mean_square(adc_input, adc_input_int, yri_cut, yri_cut_int, remainder, sim_options)

bb = 0;
vv = 0;
kk = 0;
tt = 0;
int_size_double = "single";

x3 = zeros(sim_options.Size_matrix,sim_options.Size_matrix);                       % (стр.6, (20))
x3_int = cast(zeros(sim_options.Size_matrix,sim_options.Size_matrix), int_size_double); 


DetM_2x2_multiplier_total_abs_max = cast(zeros(sim_options.num_det2x2*2,1), int_size_double);
Det2x2_sum_abs_max = cast(zeros(sim_options.num_det2x2,1), int_size_double);
Mult_DetM_3x3_array_max = cast(zeros(sim_options.num_det2x2*3,1), int_size_double);
Mult_DetM_3x3_array_mult_total_width_max = cast(zeros(sim_options.num_det2x2*3,1), int_size_double);
DetM_3x3_int_pre_sum_array_max = cast(zeros(sim_options.num_det2x2,1), int_size_double);
DetM_3x3_int_pre_sum_width_total_max = cast(zeros(sim_options.num_det2x2,1), int_size_double);
DetM_3x3_int_sum_array_max = cast(zeros(sim_options.num_det2x2,1), int_size_double);
DetM_3x3_int_sum_width_total_max = cast(zeros(sim_options.num_det2x2,1), int_size_double);
% умножители определителя 4х4
DetM_4x4_int_mult_array_max = cast(zeros(sim_options.num_det2x2*2,1), int_size_double);
% разрядность умножителей определителя 4х4
DetM_4x4_int_mult_width_total_max = cast(zeros(sim_options.num_det2x2*2,1), int_size_double);
% пресумматоры определителя 4х4
DetM_4x4_int_pre_sum_array_max = cast(zeros(sim_options.num_det2x2,1), int_size_double);
% разрядность пресумматоров определителя 4х4
DetM_4x4_int_pre_sum_width_total_max = cast(zeros(sim_options.num_det2x2,1), int_size_double);
% сумматоры определителя 4х4
DetM_4x4_int_sum_array_max = cast(zeros(5,1), int_size_double);
% разрядность сумматоров определителя 4х4
DetM_4x4_int_sum_array_width_total_max = cast(zeros(5,1), int_size_double);
% умножители определителя 5х5
DetM_5x5_int_mult_array_max = cast(zeros(5,1), int_size_double);
% разрядность умножителей 5x5
DetM_5x5_int_mult_array_width_total_max = cast(zeros(5,1), int_size_double);
% пресумматоры1 определителя 5х5
DetM_5x5_int_pre_sum1_array_max = cast(zeros(2,1), int_size_double);
% разрядность пресумматоров1 определителя 5х5
DetM_5x5_int_pre_sum1_array_width_total_max = cast(zeros(2,1), int_size_double);
% пресумматор2 определителя 5х5
DetM_5x5_int_sum3_abs_max = cast(0, int_size_double);
% разрядность пресумматора2 определителя 5х5
DetM_5x5_int_sum3_width_total_max = cast(0, int_size_double);
% сумматор определителя 5х5
DetM_5x5_int_abs_max = cast(0, int_size_double);
% разрядность сумматора определителя 5х5
DetM_5x5_int_width_total_max = cast(0, int_size_double);


%% блок для расчета первых N коэффициентов фильтра
% создаем матрицу входного сигнала
for i = 1:sim_options.Size_matrix
    x3(i,:) = adc_input(i:sim_options.Size_matrix+i-1).'; % (стр.6, (20))
	x3_int(i,:) = adc_input_int(i:sim_options.Size_matrix+i-1).'; % fi(1,12,11)
end

% рассчитываем первые N коэффициентов адаптивного фильтра
% сравнивая с задержанным сигналом ADC0 (yri_cut)
% w1 = (x3' * x3) \ x3' * yri_cut(1:N,z); % (стр.6, (19))
% w1 = lsqminnorm(x3, yri_cut(1:N,z));
   
tt = tt + 1;
%% initial determinant
det_matlab(tt) = det(x3);

[det_x3(tt), det_x3_int(tt), DetM_2x2, DetM_2x2_int, DetM_3x3_array, ...   
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
] = determinate(x3, x3_int, sim_options.int_size, sim_options.width_hilbert);

DetM_2x2_array_int(1:10) = DetM_2x2_int;
DetM_2x2_array(1:10) = DetM_2x2;
%%
DetM_3x3_array_int(1:10) = DetM_3x3_int_sum_array;
DetM_3x3_array_dd(1:10) = DetM_3x3_array;


bb = bb + 10;


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


%%
for i = 1:sim_options.Size_matrix
	kk = kk + 1;

    x3_shift = x3;
    x3_shift(1:sim_options.Size_matrix,i) = yri_cut(1:sim_options.Size_matrix);

    % int
    x3_shift_int = x3_int;
    x3_shift_int(1:sim_options.Size_matrix,i) = yri_cut_int(1:sim_options.Size_matrix); 

    %% determinant
    det_x3_shift(kk) = det(x3_shift);

    [det_out_shift(kk), det_out_shift_int(kk), DetM_2x2, DetM_2x2_int, DetM_3x3_array, ....   
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
	] = determinate(x3_shift, x3_shift_int, sim_options.int_size, sim_options.width_hilbert);

    %%
    if i == 1
       DetM_2x2_int_double = [(DetM_2x2_int(1:4))*2^-remainder;  (DetM_2x2_int(5:10))];
       DetM_3x3_int_double = [(DetM_3x3_int_sum_array(1:4)); (DetM_3x3_int_sum_array(5:10))*2^-remainder];
    elseif i == 2
       DetM_2x2_int_double = [(DetM_2x2_int(1))*2^-remainder; (DetM_2x2_int(2:4)); (DetM_2x2_int(5:7))*2^-remainder; (DetM_2x2_int(8:10))];
       DetM_3x3_int_double = [DetM_3x3_int_sum_array(1); (DetM_3x3_int_sum_array(2:4))*2^-remainder; DetM_3x3_int_sum_array(5:7); (DetM_3x3_int_sum_array(8:10))*2^-remainder];
    elseif i == 3
       DetM_2x2_int_double = [(DetM_2x2_int(1)); (DetM_2x2_int(2))*2^-remainder; (DetM_2x2_int(3:4)); (DetM_2x2_int(5))*2^-remainder; (DetM_2x2_int(6:7)); ...
       double(DetM_2x2_int(8:9))*2^-remainder; (DetM_2x2_int(10))];

       DetM_3x3_int_double = [(DetM_3x3_int_sum_array(1))*2^-remainder; (DetM_3x3_int_sum_array(2)); (DetM_3x3_int_sum_array(3:4))*2^-remainder; DetM_3x3_int_sum_array(5); 
           (DetM_3x3_int_sum_array(6:7))*2^-remainder; DetM_3x3_int_sum_array(8:9); DetM_3x3_int_sum_array(10)*2^-remainder];
    elseif i == 4
       DetM_2x2_int_double = [(DetM_2x2_int(1:2)); (DetM_2x2_int(3))*2^-remainder; (DetM_2x2_int(4:5)); (DetM_2x2_int(6))*2^-remainder; (DetM_2x2_int(7)); ...
       (DetM_2x2_int(8))*2^-remainder; (DetM_2x2_int(9)); (DetM_2x2_int(10))*2^-remainder];

       DetM_3x3_int_double = [(DetM_3x3_int_sum_array(1:2))*2^-remainder; (DetM_3x3_int_sum_array(3)); (DetM_3x3_int_sum_array(4:5))*2^-remainder; (DetM_3x3_int_sum_array(6)); (DetM_3x3_int_sum_array(7))*2^-remainder; ...
       DetM_3x3_int_sum_array(8); DetM_3x3_int_sum_array(9)*2^-remainder; DetM_3x3_int_sum_array(10)];
    elseif i == 5
       DetM_2x2_int_double = [(DetM_2x2_int(1:3)); (DetM_2x2_int(4))*2^-remainder; (DetM_2x2_int(5:6)); (DetM_2x2_int(7))*2^-remainder; (DetM_2x2_int(8)); ...
       (DetM_2x2_int(9:10))*2^-remainder];

       DetM_3x3_int_double = [(DetM_3x3_int_sum_array(1:3))*2^-remainder; (DetM_3x3_int_sum_array(4)); (DetM_3x3_int_sum_array(5:6))*2^-remainder; (DetM_3x3_int_sum_array(7)); (DetM_3x3_int_sum_array(8))*2^-remainder; ...
       (DetM_3x3_int_sum_array(9:10))];
    end

    DetM_2x2_array(bb+1:bb+10) = DetM_2x2;
    DetM_2x2_array_int(bb+1:bb+10) = DetM_2x2_int_double;
	
	DetM_3x3_array_int(bb+1:bb+10) = DetM_3x3_int_double;
	DetM_3x3_array_dd(bb+1:bb+10) = DetM_3x3_array;

    bb = bb + 10;

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
	www1(:,i) = det_x3_shift(kk) / det_matlab(tt); % double
	www1_int(:,i) = divide(det_out_shift_int(kk), det_x3_int(tt), 14); % integer

	www1_int_double(i,:) = double(www1_int(:,i))*2^-14;
end

%% filter
% умножаем входные слова на рассчитанные коэффициенты
% y_out = w1(1)*x(j+1) + w1(2)*x(j+2) + w1(3)*x(j+3) +  w1(4)*x(j+4); % Behrouz Farhang-Boroujeny, Adaptive Filters Theory and Applications  (стр. 414)

for k = 1:sim_options.Size_matrix
    dat_in_filt_double(k) = adc_input(k);
    dat_in_filt(k) = cast(adc_input_int(k),"double");
end

[y_out, y_out_int] = filter_transversal(dat_in_filt_double, www1, dat_in_filt, www1_int);

%% array out
y_array(1) = y_out;
y_array_int(1) = y_out_int;


%% part 2
for j = 1:length(yri_cut(:,1))-2*sim_options.Size_matrix
  
    % shift to left matrix input signal. Refresh matrix input signal for every new word
    for i = 1:sim_options.Size_matrix
        x3(i,:) = [x3(i,2:sim_options.Size_matrix), 0];
        x3(i,sim_options.Size_matrix) = adc_input(j+sim_options.Size_matrix-1+i); % (стр.6, (20))

        x3_int(i,:) = [x3_int(i,2:sim_options.Size_matrix), 0]; % integer
        x3_int(i,sim_options.Size_matrix) = adc_input_int(j+sim_options.Size_matrix-1+i); 
    end
    %% determinant
           
    tt = tt + 1;
    det_matlab(tt) = det(x3);

    [det_x3(tt), det_x3_int(tt), DetM_2x2, DetM_2x2_int, DetM_3x3_array,...
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
	] = determinate(x3, x3_int, sim_options.int_size, sim_options.width_hilbert); % int

    DetM_2x2_array(bb+1:bb+10) = DetM_2x2;
    DetM_2x2_array_int(bb+1:bb+10) = DetM_2x2_int;
	
	DetM_3x3_array_int(bb+1:bb+10) = DetM_3x3_int_sum_array;
	DetM_3x3_array_dd(bb+1:bb+10) = DetM_3x3_array;
	
    bb = bb + 10;
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
	
	%%
	if (det_matlab(tt) == 0)
		det_matlab(tt) = 1;
	end
            
	for i = 1:sim_options.Size_matrix
        kk = kk + 1;

		x3_shift = x3; % double
		x3_shift(1:sim_options.Size_matrix,i) = yri_cut(j+1:sim_options.Size_matrix+j); 

		x3_shift_int = x3_int; % int
		x3_shift_int(1:sim_options.Size_matrix,i) = yri_cut_int(j+1:sim_options.Size_matrix+j); 

        %%
        if (j == 625 && i == 2)
            e = 1;
        end
        %% determinant
        det_x3_shift(kk) = det(x3_shift);
        [det_out_shift(kk), det_out_shift_int(kk), DetM_2x2, DetM_2x2_int, DetM_3x3_array,...
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
        ] = determinate(x3_shift, x3_shift_int, sim_options.int_size, sim_options.width_hilbert);

        %%
        for n = 1:10
            if (DetM_2x2(n) == 21.000980680342764 )
                e = 1;
            end
        end
        %%
    if i == 1
       DetM_2x2_int_double = [(DetM_2x2_int(1:4))*2^-remainder;  (DetM_2x2_int(5:10))];
       DetM_3x3_int_double = [(DetM_3x3_int_sum_array(1:4)); (DetM_3x3_int_sum_array(5:10))*2^-remainder];
    elseif i == 2
       DetM_2x2_int_double = [(DetM_2x2_int(1))*2^-remainder; (DetM_2x2_int(2:4)); (DetM_2x2_int(5:7))*2^-remainder; (DetM_2x2_int(8:10))];
       DetM_3x3_int_double = [DetM_3x3_int_sum_array(1); (DetM_3x3_int_sum_array(2:4))*2^-remainder; DetM_3x3_int_sum_array(5:7); (DetM_3x3_int_sum_array(8:10))*2^-remainder];
    elseif i == 3
       DetM_2x2_int_double = [(DetM_2x2_int(1)); (DetM_2x2_int(2))*2^-remainder; (DetM_2x2_int(3:4)); (DetM_2x2_int(5))*2^-remainder; (DetM_2x2_int(6:7)); ...
       double(DetM_2x2_int(8:9))*2^-remainder; (DetM_2x2_int(10))];

       DetM_3x3_int_double = [(DetM_3x3_int_sum_array(1))*2^-remainder; (DetM_3x3_int_sum_array(2)); (DetM_3x3_int_sum_array(3:4))*2^-remainder; DetM_3x3_int_sum_array(5); 
           (DetM_3x3_int_sum_array(6:7))*2^-remainder; DetM_3x3_int_sum_array(8:9); DetM_3x3_int_sum_array(10)*2^-remainder];
    elseif i == 4
       DetM_2x2_int_double = [(DetM_2x2_int(1:2)); (DetM_2x2_int(3))*2^-remainder; (DetM_2x2_int(4:5)); (DetM_2x2_int(6))*2^-remainder; (DetM_2x2_int(7)); ...
       (DetM_2x2_int(8))*2^-remainder; (DetM_2x2_int(9)); (DetM_2x2_int(10))*2^-remainder];

       DetM_3x3_int_double = [(DetM_3x3_int_sum_array(1:2))*2^-remainder; (DetM_3x3_int_sum_array(3)); (DetM_3x3_int_sum_array(4:5))*2^-remainder; (DetM_3x3_int_sum_array(6)); (DetM_3x3_int_sum_array(7))*2^-remainder; ...
       DetM_3x3_int_sum_array(8); DetM_3x3_int_sum_array(9)*2^-remainder; DetM_3x3_int_sum_array(10)];
    elseif i == 5
       DetM_2x2_int_double = [(DetM_2x2_int(1:3)); (DetM_2x2_int(4))*2^-remainder; (DetM_2x2_int(5:6)); (DetM_2x2_int(7))*2^-remainder; (DetM_2x2_int(8)); ...
       (DetM_2x2_int(9:10))*2^-remainder];

       DetM_3x3_int_double = [(DetM_3x3_int_sum_array(1:3))*2^-remainder; (DetM_3x3_int_sum_array(4)); (DetM_3x3_int_sum_array(5:6))*2^-remainder; (DetM_3x3_int_sum_array(7)); (DetM_3x3_int_sum_array(8))*2^-remainder; ...
       (DetM_3x3_int_sum_array(9:10))];
    end

        DetM_2x2_array(bb+1:bb+10) = DetM_2x2;
        DetM_2x2_array_int(bb+1:bb+10) = DetM_2x2_int_double;
		
		DetM_3x3_array_int(bb+1:bb+10) = DetM_3x3_int_double;
		DetM_3x3_array_dd(bb+1:bb+10) = DetM_3x3_array;
	
        bb = bb + 10;

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

		www1(:,i) = det_x3_shift(kk) ./ det_matlab(tt); % double
			
		www1_int(:,i) = divide(det_out_shift_int(kk), det_x3_int(tt), 14); % int
		www1_int_double(i,:) = double(www1_int(:,i))*2^-14;

		%% filter
		% y_outd = 0;
		% filter input signal. Mult input words on coeff
		for k = 1:sim_options.Size_matrix
			% y_outd = y_outd + www1(k) * adc_input(j+k,z); % (стр 5, (13))
			dat_in_filt_double(k) = adc_input(j+k);
			dat_in_filt(k) = cast(adc_input_int(j+k), "double");
		end 

		[y_out, y_out_int] = filter_transversal(dat_in_filt_double, www1, dat_in_filt, www1_int);
	end

	y_array(j+1) = y_out;
	y_array_int(j+1) = y_out_int;

	end
	
    %%
    relativeError_DetM_2x2_Myfunc_double_vs_Myfunc_int = DetM_2x2_array./DetM_2x2_array_int;
    figure(8)
    plot(relativeError_DetM_2x2_Myfunc_double_vs_Myfunc_int, '-o');
    title('Относительная ошибка между определителями 2x2, найденных с помощью прямого нахождения Double vs Single')
    ylabel('Величина ошибки') 
    xlabel('Номер отсчета') 

    relativeError_DetM_3x3_Myfunc_double_vs_Myfunc_int = DetM_3x3_array_dd./DetM_3x3_array_int;
    figure(9)
    plot(relativeError_DetM_3x3_Myfunc_double_vs_Myfunc_int, '-o');
    title('Относительная ошибка между определителями 3x3, найденных с помощью прямого нахождения Double vs Single')
    ylabel('Величина ошибки') 
    xlabel('Номер отсчета') 

    %% 5x5
	det_out_shift_d = (det_out_shift_int)*2^-(remainder);
	relativeError_DetM_5x5_LU_vs_Myfunc = det_x3_shift./det_out_shift_d;
    relativeError_DetM_5x5_Myfunc_double_vs_Myfunc_int = det_out_shift./det_out_shift_d;
	figure(10)
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

    figure(11)
	subplot(2,1,1)
    plot(det_x3_shift)
    title('Значения определителей в double')
    ylabel('Значение определителя') 
    xlabel('Номер определителя') 
    subplot(2,1,2)
    plot(det_out_shift_d);
    title('Значения определителей в single')
    ylabel('Значение определителя') 
    xlabel('Номер определителя') 
end