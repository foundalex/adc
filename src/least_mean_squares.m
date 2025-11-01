function [y_array, y_array_int, DetM_2x2_array, DetM_2x2_array_int, Det2x2_mult1_abs_max, Det2x2_mult2_abs_max, Det2x2_sum_abs_max, ...
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
    ] = least_mean_square(adc_input, adc_input_int, yri_cut, yri_cut_int, M, N, width, int_size)

kk = 0;
tt = 0;
int_size_double = "double";

x3 = zeros(N,N);                       % (стр.6, (20))
x3_int = cast(zeros(N,N), int_size); 

% функция поиска определителя 2x2
num_det2x2 = 10;
Det2x2_mult1_abs_max = cast(zeros(10,1), int_size); % первые 10 умножителей1
Det2x2_mult2_abs_max = cast(zeros(10,1), int_size); % первые 10 умножителей2
Det2x2_sum_abs_max = cast(zeros(10,1), int_size); % первые 10 сумматоров

Mult_DetM_2x2_1_abs_max  = cast(zeros(3,1), int_size_double);
Mult_DetM_2x2_2_abs_max  = cast(zeros(3,1), int_size_double); 
Mult_DetM_2x2_3_abs_max  = cast(zeros(3,1), int_size_double);
Mult_DetM_2x2_4_abs_max  = cast(zeros(3,1), int_size_double);
Mult_DetM_2x2_5_abs_max  = cast(zeros(3,1), int_size_double);
Mult_DetM_2x2_6_abs_max  = cast(zeros(3,1), int_size_double);
Mult_DetM_2x2_7_abs_max  = cast(zeros(3,1), int_size_double);
Mult_DetM_2x2_8_abs_max  = cast(zeros(3,1), int_size_double);
Mult_DetM_2x2_9_abs_max  = cast(zeros(3,1), int_size_double);
Mult_DetM_2x2_10_abs_max = cast(zeros(3,1), int_size_double);   

DetM_3x3_n_11_int_sum1_abs_max = cast(0, int_size_double);
DetM_3x3_n_11_int_abs_max = cast(0, int_size_double);
DetM_3x3_n_12_int_sum1_abs_max = cast(0, int_size_double);
DetM_3x3_n_12_int_abs_max = cast(0, int_size_double);
DetM_3x3_n_13_int_sum1_abs_max = cast(0, int_size_double);
DetM_3x3_n_13_int_abs_max = cast(0, int_size_double);
DetM_3x3_n_14_int_sum1_abs_max = cast(0, int_size_double);
DetM_3x3_n_14_int_abs_max = cast(0, int_size_double);
DetM_3x3_n_22_int_sum1_abs_max = cast(0, int_size_double);
DetM_3x3_n_22_int_abs_max = cast(0, int_size_double);
DetM_3x3_n_23_int_sum1_abs_max = cast(0, int_size_double);
DetM_3x3_n_23_int_abs_max = cast(0, int_size_double);
DetM_3x3_n_24_int_sum1_abs_max = cast(0, int_size_double);
DetM_3x3_n_24_int_abs_max = cast(0, int_size_double);
DetM_3x3_n_33_int_sum1_abs_max = cast(0, int_size_double);
DetM_3x3_n_33_int_abs_max = cast(0, int_size_double);
DetM_3x3_n_34_int_sum1_abs_max = cast(0, int_size_double);
DetM_3x3_n_34_int_abs_max = cast(0, int_size_double);
DetM_3x3_n_44_int_sum1_abs_max = cast(0, int_size_double);
DetM_3x3_n_44_int_abs_max = cast(0, int_size_double);

%% блок для расчета первых N коэффициентов фильтра
%%
    
bb = 0;

% создаем матрицу входного сигнала
for i = 1:N
    x3(i,:) = adc_input(i:N+i-1).'; % (стр.6, (20))
	x3_int(i,:) = adc_input_int(i:N+i-1).'; % fi(1,12,11)
end

% рассчитываем первые N коэффициентов адаптивного фильтра
% сравнивая с задержанным сигналом ADC0 (yri_cut)
% w1 = (x3' * x3) \ x3' * yri_cut(1:N,z); % (стр.6, (19))
% w1 = lsqminnorm(x3, yri_cut(1:N,z));
   
tt = tt + 1;

%% initial determinant
det_matlab(tt) = det(x3);

[det_x3, det_x3_int, DetM_2x2, DetM_2x2_int, Det2x2_mult1_abs, Det2x2_mult2_abs, Det2x2_sum_abs, ...   
	... % первые умножители определителя 3х3
    Mult_DetM_2x2_n_1_abs, ...
    Mult_DetM_2x2_n_2_abs, ... 
    Mult_DetM_2x2_n_3_abs, ...
    Mult_DetM_2x2_n_4_abs, ...
    Mult_DetM_2x2_n_5_abs, ...
    Mult_DetM_2x2_n_6_abs, ...
    Mult_DetM_2x2_n_7_abs, ...
    Mult_DetM_2x2_n_8_abs, ...
    Mult_DetM_2x2_n_9_abs, ...
    Mult_DetM_2x2_n_10_abs, ...
	... % сумматоры определителя 3х3
	DetM_3x3_n_11_int_sum1_abs, ...
	DetM_3x3_n_11_int_abs, ...
	DetM_3x3_n_12_int_sum1_abs, ...
	DetM_3x3_n_12_int_abs, ...
	DetM_3x3_n_13_int_sum1_abs, ...
	DetM_3x3_n_13_int_abs, ...
	DetM_3x3_n_14_int_sum1_abs, ...
	DetM_3x3_n_14_int_abs, ...
	DetM_3x3_n_22_int_sum1_abs, ...
	DetM_3x3_n_22_int_abs, ...
	DetM_3x3_n_23_int_sum1_abs, ...
	DetM_3x3_n_23_int_abs, ...
	DetM_3x3_n_24_int_sum1_abs, ...
	DetM_3x3_n_24_int_abs, ...
	DetM_3x3_n_33_int_sum1_abs, ...
	DetM_3x3_n_33_int_abs, ...
	DetM_3x3_n_34_int_sum1_abs, ...
	DetM_3x3_n_34_int_abs, ...
	DetM_3x3_n_44_int_sum1_abs, ...
	DetM_3x3_n_44_int_abs ...
] = determinate(x3, x3_int, int_size, width);

%% определяем макс. значения в функции поиска определителя 2x2
for n = 1:num_det2x2
    % определяем максимальное значение на каждом из 10 умножителей
	if Det2x2_mult1_abs_max(n) < Det2x2_mult1_abs(n) 
	    Det2x2_mult1_abs_max(n) = Det2x2_mult1_abs(n);
    end

	% определяем максимальное значение на каждом из 10 умножителей
	if Det2x2_mult2_abs_max(n) < Det2x2_mult2_abs(n)
	    Det2x2_mult2_abs_max(n) = Det2x2_mult2_abs(n);
	end

	% определяем максимальное значение на каждом из 10 сумматоров
	if Det2x2_sum_abs_max(n) < Det2x2_sum_abs(n)
	    Det2x2_sum_abs_max(n) = Det2x2_sum_abs(n);
    end
end 

DetM_2x2_int_double = double(DetM_2x2_int);
relativeError_DetM_2x2 = DetM_2x2./DetM_2x2_int_double;
figure(11)
plot(relativeError_DetM_2x2, '-o');

DetM_2x2_array_int(1:10) = DetM_2x2_int_double;
DetM_2x2_array(1:10) = DetM_2x2;
bb = bb + 10;

for n = 1:3
    % определяем максимальное значение на каждом из 3 умножителей
	if Mult_DetM_2x2_1_abs_max(n) < Mult_DetM_2x2_n_1_abs(n) 
	    Mult_DetM_2x2_1_abs_max(n) = Mult_DetM_2x2_n_1_abs(n);
    end
			
	% определяем максимальное значение на каждом из 3 умножителей
	if Mult_DetM_2x2_2_abs_max(n) < Mult_DetM_2x2_n_2_abs(n) 
		Mult_DetM_2x2_2_abs_max(n) = Mult_DetM_2x2_n_2_abs(n);
	end
			
	% определяем максимальное значение на каждом из 3 умножителей
	if Mult_DetM_2x2_3_abs_max(n) < Mult_DetM_2x2_n_3_abs(n) 
		Mult_DetM_2x2_3_abs_max(n) = Mult_DetM_2x2_n_3_abs(n);
	end
			
	% определяем максимальное значение на каждом из 3 умножителей
	if Mult_DetM_2x2_4_abs_max(n) < Mult_DetM_2x2_n_4_abs(n) 
		Mult_DetM_2x2_4_abs_max(n) = Mult_DetM_2x2_n_4_abs(n);
	end
			
	% определяем максимальное значение на каждом из 3 умножителей
	if Mult_DetM_2x2_5_abs_max(n) < Mult_DetM_2x2_n_5_abs(n) 
		Mult_DetM_2x2_5_abs_max(n) = Mult_DetM_2x2_n_5_abs(n);
	end
			
	% определяем максимальное значение на каждом из 3 умножителей
	if Mult_DetM_2x2_6_abs_max(n) < Mult_DetM_2x2_n_6_abs(n) 
		Mult_DetM_2x2_6_abs_max(n) = Mult_DetM_2x2_n_6_abs(n);
	end
			
	% определяем максимальное значение на каждом из 3 умножителей
	if Mult_DetM_2x2_7_abs_max(n) < Mult_DetM_2x2_n_7_abs(n) 
		Mult_DetM_2x2_7_abs_max(n) = Mult_DetM_2x2_n_7_abs(n);
	end
			
	% определяем максимальное значение на каждом из 3 умножителей
	if Mult_DetM_2x2_8_abs_max(n) < Mult_DetM_2x2_n_8_abs(n) 
		Mult_DetM_2x2_8_abs_max(n) = Mult_DetM_2x2_n_8_abs(n);
	end
			
	% определяем максимальное значение на каждом из 3 умножителей
	if Mult_DetM_2x2_9_abs_max(n) < Mult_DetM_2x2_n_9_abs(n) 
		Mult_DetM_2x2_9_abs_max(n) = Mult_DetM_2x2_n_9_abs(n);
	end
			
	% определяем максимальное значение на каждом из 3 умножителей
	if Mult_DetM_2x2_10_abs_max(n) < Mult_DetM_2x2_n_10_abs(n) 
		Mult_DetM_2x2_10_abs_max(n) = Mult_DetM_2x2_n_10_abs(n);
	end
end

%%
% определяем максимальное значение на сумматоре определителя 3х3
if DetM_3x3_n_11_int_sum1_abs_max < DetM_3x3_n_11_int_sum1_abs 
	DetM_3x3_n_11_int_sum1_abs_max = DetM_3x3_n_11_int_sum1_abs;
end

if DetM_3x3_n_11_int_abs_max < DetM_3x3_n_11_int_abs 
	DetM_3x3_n_11_int_abs_max = DetM_3x3_n_11_int_abs;
end


if DetM_3x3_n_12_int_sum1_abs_max < DetM_3x3_n_12_int_sum1_abs 
	DetM_3x3_n_12_int_sum1_abs_max = DetM_3x3_n_12_int_sum1_abs;
end

if DetM_3x3_n_12_int_abs_max < DetM_3x3_n_12_int_abs 
	DetM_3x3_n_12_int_abs_max = DetM_3x3_n_12_int_abs;
end

if DetM_3x3_n_13_int_sum1_abs_max < DetM_3x3_n_13_int_sum1_abs 
	DetM_3x3_n_13_int_sum1_abs_max = DetM_3x3_n_13_int_sum1_abs;
end

if DetM_3x3_n_13_int_abs_max < DetM_3x3_n_13_int_abs 
	DetM_3x3_n_13_int_abs_max = DetM_3x3_n_13_int_abs;
end

if DetM_3x3_n_14_int_sum1_abs_max < DetM_3x3_n_14_int_sum1_abs 
	DetM_3x3_n_14_int_sum1_abs_max = DetM_3x3_n_14_int_sum1_abs;
end

if DetM_3x3_n_14_int_abs_max < DetM_3x3_n_14_int_abs 
	DetM_3x3_n_14_int_abs_max = DetM_3x3_n_14_int_abs;
end

if DetM_3x3_n_22_int_sum1_abs_max < DetM_3x3_n_22_int_sum1_abs 
	DetM_3x3_n_22_int_sum1_abs_max = DetM_3x3_n_22_int_sum1_abs;
end

if DetM_3x3_n_22_int_abs_max < DetM_3x3_n_22_int_abs 
	DetM_3x3_n_22_int_abs_max = DetM_3x3_n_22_int_abs;
end

if DetM_3x3_n_23_int_sum1_abs_max < DetM_3x3_n_23_int_sum1_abs 
	DetM_3x3_n_23_int_sum1_abs_max = DetM_3x3_n_23_int_sum1_abs;
end

if DetM_3x3_n_23_int_abs_max < DetM_3x3_n_23_int_abs 
	DetM_3x3_n_23_int_abs_max = DetM_3x3_n_23_int_abs;
end

if DetM_3x3_n_24_int_sum1_abs_max < DetM_3x3_n_24_int_sum1_abs 
	DetM_3x3_n_24_int_sum1_abs_max = DetM_3x3_n_24_int_sum1_abs;
end

if DetM_3x3_n_24_int_abs_max < DetM_3x3_n_24_int_abs 
	DetM_3x3_n_24_int_abs_max = DetM_3x3_n_24_int_abs;
end

if DetM_3x3_n_33_int_sum1_abs_max < DetM_3x3_n_33_int_sum1_abs 
	DetM_3x3_n_33_int_sum1_abs_max = DetM_3x3_n_33_int_sum1_abs;
end

if DetM_3x3_n_33_int_abs_max < DetM_3x3_n_33_int_abs 
	DetM_3x3_n_33_int_abs_max = DetM_3x3_n_33_int_abs;
end

if DetM_3x3_n_34_int_sum1_abs_max < DetM_3x3_n_34_int_sum1_abs 
	DetM_3x3_n_34_int_sum1_abs_max = DetM_3x3_n_34_int_sum1_abs;
end

if DetM_3x3_n_34_int_abs_max < DetM_3x3_n_34_int_abs 
	DetM_3x3_n_34_int_abs_max = DetM_3x3_n_34_int_abs;
end

if DetM_3x3_n_44_int_sum1_abs_max < DetM_3x3_n_44_int_sum1_abs 
	DetM_3x3_n_44_int_sum1_abs_max = DetM_3x3_n_44_int_sum1_abs;
end

if DetM_3x3_n_44_int_abs_max < DetM_3x3_n_44_int_abs 
	DetM_3x3_n_44_int_abs_max = DetM_3x3_n_44_int_abs;
end







%%
for i = 1:N
	kk = kk + 1;

    x3_shift = x3;
    x3_shift(1:N,i) = yri_cut(1:N);

    % int
    x3_shift_int = x3_int;
    x3_shift_int(1:N,i) = yri_cut_int(1:N); 

    %% determinant
    det_x3_shift(kk) = det(x3_shift);

    [det_out_shift(kk), det_out_shift_int(kk), DetM_2x2, DetM_2x2_int, Det2x2_mult1_abs, Det2x2_mult2_abs, Det2x2_sum_abs, ...
		Mult_DetM_2x2_n_1_abs, ...
		Mult_DetM_2x2_n_2_abs, ... 
		Mult_DetM_2x2_n_3_abs, ...
		Mult_DetM_2x2_n_4_abs, ...
		Mult_DetM_2x2_n_5_abs, ...
		Mult_DetM_2x2_n_6_abs, ...
		Mult_DetM_2x2_n_7_abs, ...
		Mult_DetM_2x2_n_8_abs, ...
		Mult_DetM_2x2_n_9_abs, ...
		Mult_DetM_2x2_n_10_abs, ...
		... %сумматоры определителя 3х3
		DetM_3x3_n_11_int_sum1_abs, ...
		DetM_3x3_n_11_int_abs, ...
		DetM_3x3_n_12_int_sum1_abs, ...
		DetM_3x3_n_12_int_abs, ...
		DetM_3x3_n_13_int_sum1_abs, ...
		DetM_3x3_n_13_int_abs, ...
		DetM_3x3_n_14_int_sum1_abs, ...
		DetM_3x3_n_14_int_abs, ...
		DetM_3x3_n_22_int_sum1_abs, ...
		DetM_3x3_n_22_int_abs, ...
		DetM_3x3_n_23_int_sum1_abs, ...
		DetM_3x3_n_23_int_abs, ...
		DetM_3x3_n_24_int_sum1_abs, ...
		DetM_3x3_n_24_int_abs, ...
		DetM_3x3_n_33_int_sum1_abs, ...
		DetM_3x3_n_33_int_abs, ...
		DetM_3x3_n_34_int_sum1_abs, ...
		DetM_3x3_n_34_int_abs, ...
		DetM_3x3_n_44_int_sum1_abs, ...
		DetM_3x3_n_44_int_abs ...
	] = determinate(x3_shift, x3_shift_int, int_size, width);
        
	%% определяем макс. значения в функции поиска определителя 2x2
	for n = 1:num_det2x2
		% определяем максимальное значение на каждом из 10 умножителей
		if Det2x2_mult1_abs_max(n) < Det2x2_mult1_abs(n) 
			Det2x2_mult1_abs_max(n) = Det2x2_mult1_abs(n);
		end

		% определяем максимальное значение на каждом из 10 умножителей
		if Det2x2_mult2_abs_max(n) < Det2x2_mult2_abs(n)
			Det2x2_mult2_abs_max(n) = Det2x2_mult2_abs(n);
		end

		% определяем максимальное значение на каждом из 10 сумматоров
		if Det2x2_sum_abs_max(n) < Det2x2_sum_abs(n)
			Det2x2_sum_abs_max(n) = Det2x2_sum_abs(n);
		end
	end
            
	for n = 1:3
		% определяем максимальное значение на каждом из 3 умножителей
		if Mult_DetM_2x2_1_abs_max(n) < Mult_DetM_2x2_n_1_abs(n) 
			Mult_DetM_2x2_1_abs_max(n) = Mult_DetM_2x2_n_1_abs(n);
		end
				
		% определяем максимальное значение на каждом из 3 умножителей
		if Mult_DetM_2x2_2_abs_max(n) < Mult_DetM_2x2_n_2_abs(n) 
			Mult_DetM_2x2_2_abs_max(n) = Mult_DetM_2x2_n_2_abs(n);
		end
				
		% определяем максимальное значение на каждом из 3 умножителей
		if Mult_DetM_2x2_3_abs_max(n) < Mult_DetM_2x2_n_3_abs(n) 
			Mult_DetM_2x2_3_abs_max(n) = Mult_DetM_2x2_n_3_abs(n);
		end
				
		% определяем максимальное значение на каждом из 3 умножителей
		if Mult_DetM_2x2_4_abs_max(n) < Mult_DetM_2x2_n_4_abs(n) 
			Mult_DetM_2x2_4_abs_max(n) = Mult_DetM_2x2_n_4_abs(n);
		end
				
		% определяем максимальное значение на каждом из 3 умножителей
		if Mult_DetM_2x2_5_abs_max(n) < Mult_DetM_2x2_n_5_abs(n) 
			Mult_DetM_2x2_5_abs_max(n) = Mult_DetM_2x2_n_5_abs(n);
		end
				
		% определяем максимальное значение на каждом из 3 умножителей
		if Mult_DetM_2x2_6_abs_max(n) < Mult_DetM_2x2_n_6_abs(n) 
			Mult_DetM_2x2_6_abs_max(n) = Mult_DetM_2x2_n_6_abs(n);
		end
				
		% определяем максимальное значение на каждом из 3 умножителей
		if Mult_DetM_2x2_7_abs_max(n) < Mult_DetM_2x2_n_7_abs(n) 
			Mult_DetM_2x2_7_abs_max(n) = Mult_DetM_2x2_n_7_abs(n);
		end
				
		% определяем максимальное значение на каждом из 3 умножителей
		if Mult_DetM_2x2_8_abs_max(n) < Mult_DetM_2x2_n_8_abs(n) 
			Mult_DetM_2x2_8_abs_max(n) = Mult_DetM_2x2_n_8_abs(n);
		end
				
		% определяем максимальное значение на каждом из 3 умножителей
		if Mult_DetM_2x2_9_abs_max(n) < Mult_DetM_2x2_n_9_abs(n) 
			Mult_DetM_2x2_9_abs_max(n) = Mult_DetM_2x2_n_9_abs(n);
		end
				
		% определяем максимальное значение на каждом из 3 умножителей
		if Mult_DetM_2x2_10_abs_max(n) < Mult_DetM_2x2_n_10_abs(n) 
			Mult_DetM_2x2_10_abs_max(n) = Mult_DetM_2x2_n_10_abs(n);
		end
	end
	
	%%
	% определяем максимальное значение на сумматоре определителя 3х3
	if DetM_3x3_n_11_int_sum1_abs_max < DetM_3x3_n_11_int_sum1_abs 
		DetM_3x3_n_11_int_sum1_abs_max = DetM_3x3_n_11_int_sum1_abs;
	end

	if DetM_3x3_n_11_int_abs_max < DetM_3x3_n_11_int_abs 
		DetM_3x3_n_11_int_abs_max = DetM_3x3_n_11_int_abs;
	end


	if DetM_3x3_n_12_int_sum1_abs_max < DetM_3x3_n_12_int_sum1_abs 
		DetM_3x3_n_12_int_sum1_abs_max = DetM_3x3_n_12_int_sum1_abs;
	end

	if DetM_3x3_n_12_int_abs_max < DetM_3x3_n_12_int_abs 
		DetM_3x3_n_12_int_abs_max = DetM_3x3_n_12_int_abs;
	end

	if DetM_3x3_n_13_int_sum1_abs_max < DetM_3x3_n_13_int_sum1_abs 
		DetM_3x3_n_13_int_sum1_abs_max = DetM_3x3_n_13_int_sum1_abs;
	end

	if DetM_3x3_n_13_int_abs_max < DetM_3x3_n_13_int_abs 
		DetM_3x3_n_13_int_abs_max = DetM_3x3_n_13_int_abs;
	end

	if DetM_3x3_n_14_int_sum1_abs_max < DetM_3x3_n_14_int_sum1_abs 
		DetM_3x3_n_14_int_sum1_abs_max = DetM_3x3_n_14_int_sum1_abs;
	end

	if DetM_3x3_n_14_int_abs_max < DetM_3x3_n_14_int_abs 
		DetM_3x3_n_14_int_abs_max = DetM_3x3_n_14_int_abs;
	end

	if DetM_3x3_n_22_int_sum1_abs_max < DetM_3x3_n_22_int_sum1_abs 
		DetM_3x3_n_22_int_sum1_abs_max = DetM_3x3_n_22_int_sum1_abs;
	end

	if DetM_3x3_n_22_int_abs_max < DetM_3x3_n_22_int_abs 
		DetM_3x3_n_22_int_abs_max = DetM_3x3_n_22_int_abs;
	end

	if DetM_3x3_n_23_int_sum1_abs_max < DetM_3x3_n_23_int_sum1_abs 
		DetM_3x3_n_23_int_sum1_abs_max = DetM_3x3_n_23_int_sum1_abs;
	end

	if DetM_3x3_n_23_int_abs_max < DetM_3x3_n_23_int_abs 
		DetM_3x3_n_23_int_abs_max = DetM_3x3_n_23_int_abs;
	end

	if DetM_3x3_n_24_int_sum1_abs_max < DetM_3x3_n_24_int_sum1_abs 
		DetM_3x3_n_24_int_sum1_abs_max = DetM_3x3_n_24_int_sum1_abs;
	end

	if DetM_3x3_n_24_int_abs_max < DetM_3x3_n_24_int_abs 
		DetM_3x3_n_24_int_abs_max = DetM_3x3_n_24_int_abs;
	end

	if DetM_3x3_n_33_int_sum1_abs_max < DetM_3x3_n_33_int_sum1_abs 
		DetM_3x3_n_33_int_sum1_abs_max = DetM_3x3_n_33_int_sum1_abs;
	end

	if DetM_3x3_n_33_int_abs_max < DetM_3x3_n_33_int_abs 
		DetM_3x3_n_33_int_abs_max = DetM_3x3_n_33_int_abs;
	end

	if DetM_3x3_n_34_int_sum1_abs_max < DetM_3x3_n_34_int_sum1_abs 
		DetM_3x3_n_34_int_sum1_abs_max = DetM_3x3_n_34_int_sum1_abs;
	end

	if DetM_3x3_n_34_int_abs_max < DetM_3x3_n_34_int_abs 
		DetM_3x3_n_34_int_abs_max = DetM_3x3_n_34_int_abs;
	end

	if DetM_3x3_n_44_int_sum1_abs_max < DetM_3x3_n_44_int_sum1_abs 
		DetM_3x3_n_44_int_sum1_abs_max = DetM_3x3_n_44_int_sum1_abs;
	end

	if DetM_3x3_n_44_int_abs_max < DetM_3x3_n_44_int_abs 
		DetM_3x3_n_44_int_abs_max = DetM_3x3_n_44_int_abs;
	end
            %%
            % if z == 3
            %     DetM_2x2_int_double = double(DetM_2x2_int)*2^-18;
            % else
            %     if i == 1
            %         DetM_2x2_int_double = [double(DetM_2x2_int(1:4))*2^-33;  double(DetM_2x2_int(5:10))];
            %     elseif i == 2
            %         DetM_2x2_int_double = [double(DetM_2x2_int(1))*2^-33; double(DetM_2x2_int(2:4)); double(DetM_2x2_int(5:7))*2^-33; double(DetM_2x2_int(8:10))];
            %     elseif i == 3
            %         DetM_2x2_int_double = [double(DetM_2x2_int(1)); double(DetM_2x2_int(2))*2^-33; double(DetM_2x2_int(3:4)); double(DetM_2x2_int(5))*2^-33; double(DetM_2x2_int(6:7)); ...
            %             double(DetM_2x2_int(8:9))*2^-33; double(DetM_2x2_int(10))];
            %     elseif i == 4
            %         DetM_2x2_int_double = [double(DetM_2x2_int(1:2)); double(DetM_2x2_int(3))*2^-33; double(DetM_2x2_int(4:5)); double(DetM_2x2_int(6))*2^-33; double(DetM_2x2_int(7)); ...
            %             double(DetM_2x2_int(8))*2^-33; double(DetM_2x2_int(9)); double(DetM_2x2_int(10))*2^-33];
            %     elseif i == 5
            %         DetM_2x2_int_double = [double(DetM_2x2_int(1:3)); double(DetM_2x2_int(4))*2^-33; double(DetM_2x2_int(5:6)); double(DetM_2x2_int(7))*2^-33; double(DetM_2x2_int(8)); ...
            %             double(DetM_2x2_int(9:10))*2^-33];
            %     end
            % end
            % 
            % relativeError_DetM_2x2 = DetM_2x2./DetM_2x2_int_double;
            % figure(11)
            % plot(relativeError_DetM_2x2, '-o');
            % title('Относительная ошибка между десятью детерминантами в double и integer')
            % xlabel('Номер детерминанта') 
            % ylabel('Величина ошибки') 
            % 
            % DetM_2x2_array_int(bb+1:bb+10) = DetM_2x2_int_double;
            % DetM_2x2_array(bb+1:bb+10) = DetM_2x2;
            % bb = bb + 10;
            %%
 

	%% divide determinant
	www1(:,i) = det_x3_shift(kk) / det_matlab(tt); % double
	www1_int(:,i) = divide(det_out_shift_int(kk), det_x3_int, 14); % integer

	www1_int_double(i,:) = double(www1_int(:,i))*2^-14;
end

%% filter
% умножаем входные слова на рассчитанные коэффициенты
% y_out = w1(1)*x(j+1) + w1(2)*x(j+2) + w1(3)*x(j+3) +  w1(4)*x(j+4); % Behrouz Farhang-Boroujeny, Adaptive Filters Theory and Applications  (стр. 414)

for k = 1:N
    dat_in_filt_double(k) = adc_input(k);
    dat_in_filt(k) = cast(adc_input_int(k),"double");
end

[y_out, y_out_int] = filter_transversal(dat_in_filt_double, www1, dat_in_filt, www1_int);

%% array out
y_array(1) = y_out;
y_array_int(1) = y_out_int;



%% part 2
for j = 1:length(yri_cut(:,1))-2*N
  
    % shift to left matrix input signal. Refresh matrix input signal for every new word
    for i = 1:N
        x3(i,:) = [x3(i,2:N), 0];
        x3(i,N) = adc_input(j+N-1+i); % (стр.6, (20))

        x3_int(i,:) = [x3_int(i,2:N), 0]; % integer
        x3_int(i,N) = adc_input_int(j+N-1+i); 
    end
    %% determinant
           
    tt = tt + 1;
    det_matlab(tt) = det(double(x3));

    [det_x3(tt), det_x3_int(tt), DetM_2x2, DetM_2x2_int, Det2x2_mult1_abs, Det2x2_mult2_abs, Det2x2_sum_abs, ...
        Mult_DetM_2x2_n_1_abs, ...
        Mult_DetM_2x2_n_2_abs, ... 
        Mult_DetM_2x2_n_3_abs, ...
        Mult_DetM_2x2_n_4_abs, ...
        Mult_DetM_2x2_n_5_abs, ...
        Mult_DetM_2x2_n_6_abs, ...
        Mult_DetM_2x2_n_7_abs, ...
        Mult_DetM_2x2_n_8_abs, ...
        Mult_DetM_2x2_n_9_abs, ...
        Mult_DetM_2x2_n_10_abs, ...    
		... % сумматоры определителя 3х3
		DetM_3x3_n_11_int_sum1_abs, ...
		DetM_3x3_n_11_int_abs, ...
		DetM_3x3_n_12_int_sum1_abs, ...
		DetM_3x3_n_12_int_abs, ...
		DetM_3x3_n_13_int_sum1_abs, ...
		DetM_3x3_n_13_int_abs, ...
		DetM_3x3_n_14_int_sum1_abs, ...
		DetM_3x3_n_14_int_abs, ...
		DetM_3x3_n_22_int_sum1_abs, ...
		DetM_3x3_n_22_int_abs, ...
		DetM_3x3_n_23_int_sum1_abs, ...
		DetM_3x3_n_23_int_abs, ...
		DetM_3x3_n_24_int_sum1_abs, ...
		DetM_3x3_n_24_int_abs, ...
		DetM_3x3_n_33_int_sum1_abs, ...
		DetM_3x3_n_33_int_abs, ...
		DetM_3x3_n_34_int_sum1_abs, ...
		DetM_3x3_n_34_int_abs, ...
		DetM_3x3_n_44_int_sum1_abs, ...
		DetM_3x3_n_44_int_abs ...
	] = determinate(x3, x3_int, int_size, width); % int

    %% определяем макс. значения в функции поиска определителя 2x2
	for n = 1:num_det2x2
		% определяем максимальное значение на каждом из 10 умножителей
		if Det2x2_mult1_abs_max(n) < Det2x2_mult1_abs(n) 
			Det2x2_mult1_abs_max(n) = Det2x2_mult1_abs(n);
		end

		% определяем максимальное значение на каждом из 10 умножителей
		if Det2x2_mult2_abs_max(n) < Det2x2_mult2_abs(n)
			Det2x2_mult2_abs_max(n) = Det2x2_mult2_abs(n);
		end

		% определяем максимальное значение на каждом из 10 сумматоров
		if Det2x2_sum_abs_max(n) < Det2x2_sum_abs(n)
			Det2x2_sum_abs_max(n) = Det2x2_sum_abs(n);
		end
	end
			
	for n = 1:3
		% определяем максимальное значение на каждом из 3 умножителей
		if Mult_DetM_2x2_1_abs_max(n) < Mult_DetM_2x2_n_1_abs(n) 
			Mult_DetM_2x2_1_abs_max(n) = Mult_DetM_2x2_n_1_abs(n);
		end
				
		% определяем максимальное значение на каждом из 3 умножителей
		if Mult_DetM_2x2_2_abs_max(n) < Mult_DetM_2x2_n_2_abs(n) 
			Mult_DetM_2x2_2_abs_max(n) = Mult_DetM_2x2_n_2_abs(n);
		end
				
		% определяем максимальное значение на каждом из 3 умножителей
		if Mult_DetM_2x2_3_abs_max(n) < Mult_DetM_2x2_n_3_abs(n) 
			Mult_DetM_2x2_3_abs_max(n) = Mult_DetM_2x2_n_3_abs(n);
		end
				
		% определяем максимальное значение на каждом из 3 умножителей
		if Mult_DetM_2x2_4_abs_max(n) < Mult_DetM_2x2_n_4_abs(n) 
			Mult_DetM_2x2_4_abs_max(n) = Mult_DetM_2x2_n_4_abs(n);
		end
				
		% определяем максимальное значение на каждом из 3 умножителей
		if Mult_DetM_2x2_5_abs_max(n) < Mult_DetM_2x2_n_5_abs(n) 
			Mult_DetM_2x2_5_abs_max(n) = Mult_DetM_2x2_n_5_abs(n);
		end
				
		% определяем максимальное значение на каждом из 3 умножителей
		if Mult_DetM_2x2_6_abs_max(n) < Mult_DetM_2x2_n_6_abs(n) 
			Mult_DetM_2x2_6_abs_max(n) = Mult_DetM_2x2_n_6_abs(n);
		end
				
		% определяем максимальное значение на каждом из 3 умножителей
		if Mult_DetM_2x2_7_abs_max(n) < Mult_DetM_2x2_n_7_abs(n) 
			Mult_DetM_2x2_7_abs_max(n) = Mult_DetM_2x2_n_7_abs(n);
		end
				
		% определяем максимальное значение на каждом из 3 умножителей
		if Mult_DetM_2x2_8_abs_max(n) < Mult_DetM_2x2_n_8_abs(n) 
			Mult_DetM_2x2_8_abs_max(n) = Mult_DetM_2x2_n_8_abs(n);
		end
				
		% определяем максимальное значение на каждом из 3 умножителей
		if Mult_DetM_2x2_9_abs_max(n) < Mult_DetM_2x2_n_9_abs(n) 
			Mult_DetM_2x2_9_abs_max(n) = Mult_DetM_2x2_n_9_abs(n);
		end
				
		% определяем максимальное значение на каждом из 3 умножителей
		if Mult_DetM_2x2_10_abs_max(n) < Mult_DetM_2x2_n_10_abs(n) 
			Mult_DetM_2x2_10_abs_max(n) = Mult_DetM_2x2_n_10_abs(n);
		end
	end
	
	%%
	% определяем максимальное значение на сумматоре определителя 3х3
	if DetM_3x3_n_11_int_sum1_abs_max < DetM_3x3_n_11_int_sum1_abs 
		DetM_3x3_n_11_int_sum1_abs_max = DetM_3x3_n_11_int_sum1_abs;
	end
	
	if DetM_3x3_n_11_int_abs_max < DetM_3x3_n_11_int_abs 
		DetM_3x3_n_11_int_abs_max = DetM_3x3_n_11_int_abs;
	end
	
	
	if DetM_3x3_n_12_int_sum1_abs_max < DetM_3x3_n_12_int_sum1_abs 
		DetM_3x3_n_12_int_sum1_abs_max = DetM_3x3_n_12_int_sum1_abs;
	end
	
	if DetM_3x3_n_12_int_abs_max < DetM_3x3_n_12_int_abs 
		DetM_3x3_n_12_int_abs_max = DetM_3x3_n_12_int_abs;
	end
	
	if DetM_3x3_n_13_int_sum1_abs_max < DetM_3x3_n_13_int_sum1_abs 
		DetM_3x3_n_13_int_sum1_abs_max = DetM_3x3_n_13_int_sum1_abs;
	end
	
	if DetM_3x3_n_13_int_abs_max < DetM_3x3_n_13_int_abs 
		DetM_3x3_n_13_int_abs_max = DetM_3x3_n_13_int_abs;
	end
	
	if DetM_3x3_n_14_int_sum1_abs_max < DetM_3x3_n_14_int_sum1_abs 
		DetM_3x3_n_14_int_sum1_abs_max = DetM_3x3_n_14_int_sum1_abs;
	end
	
	if DetM_3x3_n_14_int_abs_max < DetM_3x3_n_14_int_abs 
		DetM_3x3_n_14_int_abs_max = DetM_3x3_n_14_int_abs;
	end
	
	if DetM_3x3_n_22_int_sum1_abs_max < DetM_3x3_n_22_int_sum1_abs 
		DetM_3x3_n_22_int_sum1_abs_max = DetM_3x3_n_22_int_sum1_abs;
	end
	
	if DetM_3x3_n_22_int_abs_max < DetM_3x3_n_22_int_abs 
		DetM_3x3_n_22_int_abs_max = DetM_3x3_n_22_int_abs;
	end
	
	if DetM_3x3_n_23_int_sum1_abs_max < DetM_3x3_n_23_int_sum1_abs 
		DetM_3x3_n_23_int_sum1_abs_max = DetM_3x3_n_23_int_sum1_abs;
	end
	
	if DetM_3x3_n_23_int_abs_max < DetM_3x3_n_23_int_abs 
		DetM_3x3_n_23_int_abs_max = DetM_3x3_n_23_int_abs;
	end
	
	if DetM_3x3_n_24_int_sum1_abs_max < DetM_3x3_n_24_int_sum1_abs 
		DetM_3x3_n_24_int_sum1_abs_max = DetM_3x3_n_24_int_sum1_abs;
	end
	
	if DetM_3x3_n_24_int_abs_max < DetM_3x3_n_24_int_abs 
		DetM_3x3_n_24_int_abs_max = DetM_3x3_n_24_int_abs;
	end
	
	if DetM_3x3_n_33_int_sum1_abs_max < DetM_3x3_n_33_int_sum1_abs 
		DetM_3x3_n_33_int_sum1_abs_max = DetM_3x3_n_33_int_sum1_abs;
	end
	
	if DetM_3x3_n_33_int_abs_max < DetM_3x3_n_33_int_abs 
		DetM_3x3_n_33_int_abs_max = DetM_3x3_n_33_int_abs;
	end
	
	if DetM_3x3_n_34_int_sum1_abs_max < DetM_3x3_n_34_int_sum1_abs 
		DetM_3x3_n_34_int_sum1_abs_max = DetM_3x3_n_34_int_sum1_abs;
	end
	
	if DetM_3x3_n_34_int_abs_max < DetM_3x3_n_34_int_abs 
		DetM_3x3_n_34_int_abs_max = DetM_3x3_n_34_int_abs;
	end
	
	if DetM_3x3_n_44_int_sum1_abs_max < DetM_3x3_n_44_int_sum1_abs 
		DetM_3x3_n_44_int_sum1_abs_max = DetM_3x3_n_44_int_sum1_abs;
	end
	
	if DetM_3x3_n_44_int_abs_max < DetM_3x3_n_44_int_abs 
		DetM_3x3_n_44_int_abs_max = DetM_3x3_n_44_int_abs;
	end

	%%
	
	DetM_2x2_int_double = double(DetM_2x2_int);
	% relativeError_DetM_2x2 = DetM_2x2./DetM_2x2_int_double;
	% figure(10)
	% plot(relativeError_DetM_2x2, '-o');

	DetM_2x2_array_int(bb+1:bb+10) = DetM_2x2_int_double;
	DetM_2x2_array(bb+1:bb+10) = DetM_2x2;
	bb = bb + 10;
	%%
	if (det_matlab(tt) == 0)
		det_matlab(tt) = 1;
	end
            
	for i = 1:N
        kk = kk + 1;

		x3_shift = x3; % double
		x3_shift(1:N,i) = yri_cut(j+1:N+j); 

		x3_shift_int = x3_int; % int
		x3_shift_int(1:N,i) = yri_cut_int(j+1:N+j); 

        %% determinant
        det_x3_shift(kk) = det(x3_shift);
        [det_out_shift(kk), det_out_shift_int(kk), DetM_2x2, DetM_2x2_int, Det2x2_mult1_abs, Det2x2_mult2_abs, Det2x2_sum_abs, ...
			Mult_DetM_2x2_n_1_abs, ...
			Mult_DetM_2x2_n_2_abs, ... 
			Mult_DetM_2x2_n_3_abs, ...
			Mult_DetM_2x2_n_4_abs, ...
			Mult_DetM_2x2_n_5_abs, ...
			Mult_DetM_2x2_n_6_abs, ...
			Mult_DetM_2x2_n_7_abs, ...
			Mult_DetM_2x2_n_8_abs, ...
			Mult_DetM_2x2_n_9_abs, ...
			Mult_DetM_2x2_n_10_abs, ...    
			... % сумматоры определителя 3х3
			DetM_3x3_n_11_int_sum1_abs, ...
			DetM_3x3_n_11_int_abs, ...
			DetM_3x3_n_12_int_sum1_abs, ...
			DetM_3x3_n_12_int_abs, ...
			DetM_3x3_n_13_int_sum1_abs, ...
			DetM_3x3_n_13_int_abs, ...
			DetM_3x3_n_14_int_sum1_abs, ...
			DetM_3x3_n_14_int_abs, ...
			DetM_3x3_n_22_int_sum1_abs, ...
			DetM_3x3_n_22_int_abs, ...
			DetM_3x3_n_23_int_sum1_abs, ...
			DetM_3x3_n_23_int_abs, ...
			DetM_3x3_n_24_int_sum1_abs, ...
			DetM_3x3_n_24_int_abs, ...
			DetM_3x3_n_33_int_sum1_abs, ...
			DetM_3x3_n_33_int_abs, ...
			DetM_3x3_n_34_int_sum1_abs, ...
			DetM_3x3_n_34_int_abs, ...
			DetM_3x3_n_44_int_sum1_abs, ...
			DetM_3x3_n_44_int_abs ...
        ] = determinate(x3_shift, x3_shift_int, int_size, width);

		%% определяем макс. значения в функции поиска определителя 2x2
		for n = 1:num_det2x2
			% определяем максимальное значение на каждом из 10 умножителей
			if Det2x2_mult1_abs_max(n) < Det2x2_mult1_abs(n) 
				Det2x2_mult1_abs_max(n) = Det2x2_mult1_abs(n);
			end

			% определяем максимальное значение на каждом из 10 умножителей
			if Det2x2_mult2_abs_max(n) < Det2x2_mult2_abs(n)
				Det2x2_mult2_abs_max(n) = Det2x2_mult2_abs(n);
			end

			% определяем максимальное значение на каждом из 10 сумматоров
			if Det2x2_sum_abs_max(n) < Det2x2_sum_abs(n)
				Det2x2_sum_abs_max(n) = Det2x2_sum_abs(n);
			end
        end
				
		for n = 1:3
			% определяем максимальное значение на каждом из 3 умножителей
			if Mult_DetM_2x2_1_abs_max(n) < Mult_DetM_2x2_n_1_abs(n) 
				Mult_DetM_2x2_1_abs_max(n) = Mult_DetM_2x2_n_1_abs(n);
			end
					
			% определяем максимальное значение на каждом из 3 умножителей
			if Mult_DetM_2x2_2_abs_max(n) < Mult_DetM_2x2_n_2_abs(n) 
				Mult_DetM_2x2_2_abs_max(n) = Mult_DetM_2x2_n_2_abs(n);
			end
					
			% определяем максимальное значение на каждом из 3 умножителей
			if Mult_DetM_2x2_3_abs_max(n) < Mult_DetM_2x2_n_3_abs(n) 
				Mult_DetM_2x2_3_abs_max(n) = Mult_DetM_2x2_n_3_abs(n);
			end
					
			% определяем максимальное значение на каждом из 3 умножителей
			if Mult_DetM_2x2_4_abs_max(n) < Mult_DetM_2x2_n_4_abs(n) 
				Mult_DetM_2x2_4_abs_max(n) = Mult_DetM_2x2_n_4_abs(n);
			end
					
			% определяем максимальное значение на каждом из 3 умножителей
			if Mult_DetM_2x2_5_abs_max(n) < Mult_DetM_2x2_n_5_abs(n) 
				Mult_DetM_2x2_5_abs_max(n) = Mult_DetM_2x2_n_5_abs(n);
			end
					
			% определяем максимальное значение на каждом из 3 умножителей
			if Mult_DetM_2x2_6_abs_max(n) < Mult_DetM_2x2_n_6_abs(n) 
				Mult_DetM_2x2_6_abs_max(n) = Mult_DetM_2x2_n_6_abs(n);
			end
					
			% определяем максимальное значение на каждом из 3 умножителей
			if Mult_DetM_2x2_7_abs_max(n) < Mult_DetM_2x2_n_7_abs(n) 
				Mult_DetM_2x2_7_abs_max(n) = Mult_DetM_2x2_n_7_abs(n);
			end
					
			% определяем максимальное значение на каждом из 3 умножителей
			if Mult_DetM_2x2_8_abs_max(n) < Mult_DetM_2x2_n_8_abs(n) 
				Mult_DetM_2x2_8_abs_max(n) = Mult_DetM_2x2_n_8_abs(n);
			end
					
			% определяем максимальное значение на каждом из 3 умножителей
			if Mult_DetM_2x2_9_abs_max(n) < Mult_DetM_2x2_n_9_abs(n) 
				Mult_DetM_2x2_9_abs_max(n) = Mult_DetM_2x2_n_9_abs(n);
			end
					
			% определяем максимальное значение на каждом из 3 умножителей
			if Mult_DetM_2x2_10_abs_max(n) < Mult_DetM_2x2_n_10_abs(n) 
				Mult_DetM_2x2_10_abs_max(n) = Mult_DetM_2x2_n_10_abs(n);
			end
		end
		
		%%
		% определяем максимальное значение на сумматоре определителя 3х3
		if DetM_3x3_n_11_int_sum1_abs_max < DetM_3x3_n_11_int_sum1_abs 
			DetM_3x3_n_11_int_sum1_abs_max = DetM_3x3_n_11_int_sum1_abs;
		end

		if DetM_3x3_n_11_int_abs_max < DetM_3x3_n_11_int_abs 
			DetM_3x3_n_11_int_abs_max = DetM_3x3_n_11_int_abs;
		end
		
		
		if DetM_3x3_n_12_int_sum1_abs_max < DetM_3x3_n_12_int_sum1_abs 
			DetM_3x3_n_12_int_sum1_abs_max = DetM_3x3_n_12_int_sum1_abs;
		end
		
		if DetM_3x3_n_12_int_abs_max < DetM_3x3_n_12_int_abs 
			DetM_3x3_n_12_int_abs_max = DetM_3x3_n_12_int_abs;
		end
		
		if DetM_3x3_n_13_int_sum1_abs_max < DetM_3x3_n_13_int_sum1_abs 
			DetM_3x3_n_13_int_sum1_abs_max = DetM_3x3_n_13_int_sum1_abs;
		end
		
		if DetM_3x3_n_13_int_abs_max < DetM_3x3_n_13_int_abs 
			DetM_3x3_n_13_int_abs_max = DetM_3x3_n_13_int_abs;
		end
		
		if DetM_3x3_n_14_int_sum1_abs_max < DetM_3x3_n_14_int_sum1_abs 
			DetM_3x3_n_14_int_sum1_abs_max = DetM_3x3_n_14_int_sum1_abs;
		end
		
		if DetM_3x3_n_14_int_abs_max < DetM_3x3_n_14_int_abs 
			DetM_3x3_n_14_int_abs_max = DetM_3x3_n_14_int_abs;
		end
		
		if DetM_3x3_n_22_int_sum1_abs_max < DetM_3x3_n_22_int_sum1_abs 
			DetM_3x3_n_22_int_sum1_abs_max = DetM_3x3_n_22_int_sum1_abs;
		end
		
		if DetM_3x3_n_22_int_abs_max < DetM_3x3_n_22_int_abs 
			DetM_3x3_n_22_int_abs_max = DetM_3x3_n_22_int_abs;
		end
		
		if DetM_3x3_n_23_int_sum1_abs_max < DetM_3x3_n_23_int_sum1_abs 
			DetM_3x3_n_23_int_sum1_abs_max = DetM_3x3_n_23_int_sum1_abs;
		end
		
		if DetM_3x3_n_23_int_abs_max < DetM_3x3_n_23_int_abs 
			DetM_3x3_n_23_int_abs_max = DetM_3x3_n_23_int_abs;
		end
		
		if DetM_3x3_n_24_int_sum1_abs_max < DetM_3x3_n_24_int_sum1_abs 
			DetM_3x3_n_24_int_sum1_abs_max = DetM_3x3_n_24_int_sum1_abs;
		end
		
		if DetM_3x3_n_24_int_abs_max < DetM_3x3_n_24_int_abs 
			DetM_3x3_n_24_int_abs_max = DetM_3x3_n_24_int_abs;
		end
		
		if DetM_3x3_n_33_int_sum1_abs_max < DetM_3x3_n_33_int_sum1_abs 
			DetM_3x3_n_33_int_sum1_abs_max = DetM_3x3_n_33_int_sum1_abs;
		end
		
		if DetM_3x3_n_33_int_abs_max < DetM_3x3_n_33_int_abs 
			DetM_3x3_n_33_int_abs_max = DetM_3x3_n_33_int_abs;
		end
		
		if DetM_3x3_n_34_int_sum1_abs_max < DetM_3x3_n_34_int_sum1_abs 
			DetM_3x3_n_34_int_sum1_abs_max = DetM_3x3_n_34_int_sum1_abs;
		end
		
		if DetM_3x3_n_34_int_abs_max < DetM_3x3_n_34_int_abs 
			DetM_3x3_n_34_int_abs_max = DetM_3x3_n_34_int_abs;
		end
		
		if DetM_3x3_n_44_int_sum1_abs_max < DetM_3x3_n_44_int_sum1_abs 
			DetM_3x3_n_44_int_sum1_abs_max = DetM_3x3_n_44_int_sum1_abs;
		end
		
		if DetM_3x3_n_44_int_abs_max < DetM_3x3_n_44_int_abs 
			DetM_3x3_n_44_int_abs_max = DetM_3x3_n_44_int_abs;
		end
				%%

                % if z == 3
                %     if i == 1
                %         DetM_2x2_int_double = [double(DetM_2x2_int(1:4))*2^-18;  double(DetM_2x2_int(5:10))];
                %     elseif i == 2
                %         DetM_2x2_int_double = [double(DetM_2x2_int(1))*2^-18; double(DetM_2x2_int(2:4)); double(DetM_2x2_int(5:7))*2^-18; double(DetM_2x2_int(8:10))];
                %     elseif i == 3
                %         DetM_2x2_int_double = [double(DetM_2x2_int(1)); double(DetM_2x2_int(2))*2^-18; double(DetM_2x2_int(3:4)); double(DetM_2x2_int(5))*2^-18; double(DetM_2x2_int(6:7)); ...
                %         double(DetM_2x2_int(8:9))*2^-18; double(DetM_2x2_int(10))];
                %     elseif i == 4
                %         DetM_2x2_int_double = [double(DetM_2x2_int(1:2)); double(DetM_2x2_int(3))*2^-18; double(DetM_2x2_int(4:5)); double(DetM_2x2_int(6))*2^-18; double(DetM_2x2_int(7)); ...
                %         double(DetM_2x2_int(8))*2^-18; double(DetM_2x2_int(9)); double(DetM_2x2_int(10))*2^-18];
                %     elseif i == 5
                %         DetM_2x2_int_double = [double(DetM_2x2_int(1:3)); double(DetM_2x2_int(4))*2^-18; double(DetM_2x2_int(5:6)); double(DetM_2x2_int(7))*2^-18; double(DetM_2x2_int(8)); ...
                %         double(DetM_2x2_int(9:10))*2^-18];
                %     end
                % else
                %     if i == 1
                %         DetM_2x2_int_double = [double(DetM_2x2_int(1:4))*2^-30;  double(DetM_2x2_int(5:10))];
                %     elseif i == 2
                %         DetM_2x2_int_double = [double(DetM_2x2_int(1))*2^-30; double(DetM_2x2_int(2:4)); double(DetM_2x2_int(5:7))*2^-30; double(DetM_2x2_int(8:10))];
                %     elseif i == 3
                %         DetM_2x2_int_double = [double(DetM_2x2_int(1)); double(DetM_2x2_int(2))*2^-30; double(DetM_2x2_int(3:4)); double(DetM_2x2_int(5))*2^-30; double(DetM_2x2_int(6:7)); ...
                %         double(DetM_2x2_int(8:9))*2^-30; double(DetM_2x2_int(10))];
                %     elseif i == 4
                %         DetM_2x2_int_double = [double(DetM_2x2_int(1:2)); double(DetM_2x2_int(3))*2^-30; double(DetM_2x2_int(4:5)); double(DetM_2x2_int(6))*2^-30; double(DetM_2x2_int(7)); ...
                %         double(DetM_2x2_int(8))*2^-30; double(DetM_2x2_int(9)); double(DetM_2x2_int(10))*2^-30];
                %     elseif i == 5
                %         DetM_2x2_int_double = [double(DetM_2x2_int(1:3)); double(DetM_2x2_int(4))*2^-30; double(DetM_2x2_int(5:6)); double(DetM_2x2_int(7))*2^-30; double(DetM_2x2_int(8)); ...
                %         double(DetM_2x2_int(9:10))*2^-30];
                %     end
                % end

                % relativeError_DetM_2x2 = DetM_2x2./DetM_2x2_int_double;
                % figure(10)
                % plot(relativeError_DetM_2x2, '-o');

                % DetM_2x2_array_int(bb+1:bb+10,z-1) = DetM_2x2_int_double;
                % DetM_2x2_array(bb+1:bb+10,z-1) = DetM_2x2;
                % bb = bb + 10;
                %%

                %% divide determinant

		www1(:,i) = det_x3_shift(kk) ./ det_matlab(tt); % double
	
		www1_int(:,i) = divide(det_out_shift_int(kk), det_x3_int(tt), 14); % int
		www1_int_double(i,:) = double(www1_int(:,i))*2^-14;

		%% filter
		% y_outd = 0;
		% filter input signal. Mult input words on coeff
		for k = 1:N
			% y_outd = y_outd + www1(k) * adc_input(j+k,z); % (стр 5, (13))
			dat_in_filt_double(k) = adc_input(j+k);
			dat_in_filt(k) = cast(adc_input_int(j+k), "double");
		end 

		[y_out, y_out_int] = filter_transversal(dat_in_filt_double, www1, dat_in_filt, www1_int);
	end

	y_array(j+1) = y_out;
	y_array_int(j+1) = y_out_int;

end
end