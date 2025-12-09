function [y_array, y_array_double, y_array_int, determinate_struct, Divide_max, adaptive_filter_struct_max ...
    ] = least_mean_square(adc_input, yri_cut, yri_cut_int, read_max_width, sim_options)

bb = 0;
vv = 0;
kk = 0;
nn = 0;

x3 = cast(zeros(sim_options.Size_matrix,sim_options.Size_matrix), "double");                     
x3_int = cast(zeros(sim_options.Size_matrix,sim_options.Size_matrix), sim_options.type_fir_out); 

determinate_struct = struct;
determinate_struct.DetM_2x2_multiplier_total_abs_max = cast(zeros(sim_options.num_det2x2*2,1), sim_options.type_2x2_det);
determinate_struct.Det2x2_sum_abs_max = cast(zeros(sim_options.num_det2x2,1), sim_options.type_2x2_det);

determinate_struct.Mult_DetM_3x3_array_max = cast(zeros(sim_options.num_det2x2*3,1), sim_options.type_3x3_det);
determinate_struct.Mult_DetM_3x3_array_mult_total_width_max = cast(zeros(sim_options.num_det2x2*3,1), sim_options.type_3x3_det);
determinate_struct.DetM_3x3_int_pre_sum_array_max = cast(zeros(sim_options.num_det2x2,1), sim_options.type_3x3_det);
determinate_struct.DetM_3x3_int_pre_sum_width_total_max = cast(zeros(sim_options.num_det2x2,1), sim_options.type_3x3_det);
determinate_struct.DetM_3x3_int_sum_array_max = cast(zeros(sim_options.num_det2x2,1), sim_options.type_3x3_det);
determinate_struct.DetM_3x3_int_sum_width_total_max = cast(zeros(sim_options.num_det2x2,1), sim_options.type_3x3_det);
% умножители определителя 4х4
determinate_struct.DetM_4x4_int_mult_array_max = cast(zeros(sim_options.num_det2x2*2,1), sim_options.type_4x4_det);
% разрядность умножителей определителя 4х4
determinate_struct.DetM_4x4_int_mult_width_total_max = cast(zeros(sim_options.num_det2x2*2,1), sim_options.type_4x4_det);
% пресумматоры определителя 4х4
determinate_struct.DetM_4x4_int_pre_sum_array_max = cast(zeros(sim_options.num_det2x2,1), sim_options.type_4x4_det);
% разрядность пресумматоров определителя 4х4
determinate_struct.DetM_4x4_int_pre_sum_width_total_max = cast(zeros(sim_options.num_det2x2,1), sim_options.type_4x4_det);
% сумматоры определителя 4х4
determinate_struct.DetM_4x4_int_sum_array_max = cast(zeros(sim_options.Size_matrix,1), sim_options.type_4x4_det);
% разрядность сумматоров определителя 4х4
determinate_struct.DetM_4x4_int_sum_array_width_total_max = cast(zeros(sim_options.Size_matrix,1), sim_options.type_4x4_det);
% умножители определителя 5х5
determinate_struct.DetM_5x5_int_mult_array_max = cast(zeros(sim_options.Size_matrix,1), sim_options.type_5x5_det);
% разрядность умножителей 5x5
determinate_struct.DetM_5x5_int_mult_array_width_total_max = cast(zeros(sim_options.Size_matrix,1), sim_options.type_5x5_det);
% пресумматоры1 определителя 5х5
determinate_struct.DetM_5x5_int_pre_sum1_array_max = cast(zeros(2,1), sim_options.type_5x5_det);
% разрядность пресумматоров1 определителя 5х5
determinate_struct.DetM_5x5_int_pre_sum1_array_width_total_max = cast(zeros(2,1), sim_options.type_5x5_det);
% пресумматор2 определителя 5х5
determinate_struct.DetM_5x5_int_sum3_abs_max = cast(0, sim_options.type_5x5_det);
% разрядность пресумматора2 определителя 5х5
determinate_struct.DetM_5x5_int_sum3_width_total_max = cast(0, sim_options.type_5x5_det);
% сумматор определителя 5х5
determinate_struct.DetM_5x5_int_abs_max = cast(0, sim_options.type_5x5_det);
% разрядность сумматора определителя 5х5
determinate_struct.DetM_5x5_int_width_total_max = cast(0, sim_options.type_5x5_det);
% начальный определитель
determinate_struct.Det_x3_int_max = cast(0, sim_options.type_5x5_det);
% выход делителя
Divide_max = cast(0, sim_options.type_divide_out);

www1 = zeros(sim_options.Size_matrix,1);
www1_double = zeros(sim_options.Size_matrix,1);
www1_double_abs = zeros(sim_options.Size_matrix,1);
www1_int = cast(zeros(sim_options.Size_matrix,1), sim_options.type_divide_out);
www1_int_abs = cast(zeros(sim_options.Size_matrix,1), sim_options.type_divide_out);

% Структура адаптивного фильтра
adaptive_filter_struct_max = struct;
adaptive_filter_struct_max.Adaptive_filter_mult_array_max = cast(zeros(sim_options.Size_matrix*sim_options.Size_matrix,1),sim_options.type_mult_in_adaptive_filter);
adaptive_filter_struct_max.Adaptive_filter_mult_total_width = cast(zeros(sim_options.Size_matrix*sim_options.Size_matrix,1),sim_options.type_mult_in_adaptive_filter);
adaptive_filter_struct_max.Adaptive_filter_sum_array_max = cast(zeros(sim_options.Size_matrix*2,3),sim_options.type_add_in_adaptive_filter);
adaptive_filter_struct_max.Adaptive_filter_sum_total_width = cast(zeros(sim_options.Size_matrix*2,3),sim_options.type_add_in_adaptive_filter);

for j = 1:sim_options.Size_matrix:length(yri_cut(:,1))-sim_options.Size_matrix+1 

    % Создаем матрицу 5х5 из отсчетов сигнала
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

    %% Поиск определителя матрицы, состоящей из отсчетов сигнала с i-го суб-АЦП
    det_matlab(j) = det(x3);
    if (det_matlab(j) == 0)
        det_matlab(j) = 1;
        disp('Determinant LU x3 equal 0');
    end

    % x31 = x3 * x3';
	[det_x3(j), det_x3_int(j), DetM_2x2, Det_2x2_LU_matlab, DetM_3x3_array, Det_3x3_LU_matlab, DetM_4x4_array, Det_4x4_LU_matlab, s...
        ] = determinate(x3, x3_int, read_max_width.det_max_width, sim_options); % int

    % det_x3_int(j) = bitshift(det_x3_int(j),-2);
    % det_x31 = sqrt(det_x3(j));

	DetM_2x2_array(bb+1:bb+10) = DetM_2x2;
	DetM_2x2_array_int(bb+1:bb+10) = s.Det2x2_sum_abs;
	Det_2x2_LU_matlab_array(bb+1:bb+10) = Det_2x2_LU_matlab;
	
	DetM_3x3_array_int(bb+1:bb+10) = s.DetM_3x3_int_sum_array;
	DetM_3x3_array_dd(bb+1:bb+10) = DetM_3x3_array;
	Det_3x3_LU_matlab_array(bb+1:bb+10) = Det_3x3_LU_matlab;

	DetM_4x4_array_int(vv+1:vv+5) = s.DetM_4x4_int_sum_array;
	DetM_4x4_array_dd(vv+1:vv+5) = DetM_4x4_array;
	Det_4x4_LU_matlab_array(vv+1:vv+5) = Det_4x4_LU_matlab;

	bb = bb + 10;
	vv = vv + 5;

	% Запись в структуру максимальных значений сумматоров и множителей
    % определителя
	determinate_struct = compare_determinante(s, determinate_struct, sim_options);
	
	for i = 1:sim_options.Size_matrix
        % Добавляем в матрицу столбец эталонного сигнала с фильтра дробной
        % задержки
		kk = kk + 1;

		x3_shift = x3;
		x3_shift(1:sim_options.Size_matrix,i) = yri_cut(j:j+sim_options.Size_matrix-1);

        % x3_shift1 = x3_shift * x3_shift';

		x3_shift_int = x3_int;
		x3_shift_int(1:sim_options.Size_matrix,i) = yri_cut_int(j:j+sim_options.Size_matrix-1); 

		w1 = lsqminnorm(x3, yri_cut(j:sim_options.Size_matrix+j-1));
		%% Поиск определителя

        % функция матлаб
		det_x3_shift(kk) = det(x3_shift);
        if (det_x3_shift(kk) == 0)
            det_x3_shift(kk) = 1;
            disp('Determinant LU x3_shift equal 0');
        end

        % собственная функция
		[det_out_shift(kk), det_out_shift_int(kk), DetM_2x2, Det_2x2_LU_matlab, DetM_3x3_array, Det_3x3_LU_matlab, DetM_4x4_array, Det_4x4_LU_matlab, s...
            ] = determinate(x3_shift, x3_shift_int, read_max_width.det_max_width, sim_options);

        % det_out_shift_int(kk) = bitshift(det_out_shift_int(kk),-2); 
        % det_x31_shift = sqrt(det_out_shift(kk));

		DetM_2x2_array(bb+1:bb+10) = DetM_2x2;
		DetM_2x2_array_int(bb+1:bb+10) = s.Det2x2_sum_abs;
		Det_2x2_LU_matlab_array(bb+1:bb+10) = Det_2x2_LU_matlab;
	
		DetM_3x3_array_int(bb+1:bb+10) = s.DetM_3x3_int_sum_array;
		DetM_3x3_array_dd(bb+1:bb+10) = DetM_3x3_array;
		Det_3x3_LU_matlab_array(bb+1:bb+10) = Det_3x3_LU_matlab;
	
		DetM_4x4_array_int(vv+1:vv+5) = s.DetM_4x4_int_sum_array;
		DetM_4x4_array_dd(vv+1:vv+5) = DetM_4x4_array;
		Det_4x4_LU_matlab_array(vv+1:vv+5) = Det_4x4_LU_matlab;

		bb = bb + 10;
		vv = vv + 5;

        % Запись в структуру максимальных значений сумматоров и множителей
        % определителя
		determinate_struct = compare_determinante(s, determinate_struct, sim_options);

		%% Деление определителя матрицы с эталонным сигналом на определитель матрицы сигнала с i-го суб-АЦП
        % для получения коэффициентов адаптивного фильтра
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

	%% Адаптивный фильтр
    y_outd = w1(1).*x3(:,1)+w1(2).*x3(:,2)+w1(3).*x3(:,3)+w1(4).*x3(:,4)+w1(5).*x3(:,5);

    x3_int_c = cast(x3_int, sim_options.type_mult_in_adaptive_filter);
    www1_int_c = cast(www1_int, sim_options.type_mult_in_adaptive_filter);

    [y_out, y_out_int, adaptive_filter_structure] = adaptive_filter(x3, www1_double, x3_int_c, www1_int_c, read_max_width.adaptive_max_width, sim_options);

    % data_outd = www1_double(1).*x31(:,1)+www1_double(2).*x31(:,2)+www1_double(3).*x31(:,3)+www1_double(4).*x31(:,4)+www1_double(5).*x31(:,5);
    y_out_double =round(y_out * 2^-sim_options.divide_factor);

    % округление значений после фильтра
    y_out_int_shift = round_int(y_out_int, sim_options.divide_factor, sim_options.type_fir_out);


	%% определяем макс. значения
	for n = 1:sim_options.Size_matrix*sim_options.Size_matrix
		% определяем максимальное значение на каждом умножителе
		if adaptive_filter_struct_max.Adaptive_filter_mult_array_max(n) < adaptive_filter_structure.mult_int_abs(n) 
			adaptive_filter_struct_max.Adaptive_filter_mult_array_max(n) = adaptive_filter_structure.mult_int_abs(n); 
        end
		% определяем максимальную разрядность умножителей
		if adaptive_filter_struct_max.Adaptive_filter_mult_total_width(n) < adaptive_filter_structure.mult_int_total_width(n) 
			adaptive_filter_struct_max.Adaptive_filter_mult_total_width(n) = adaptive_filter_structure.mult_int_total_width(n);
        end
    end

    for n = 1:sim_options.Size_matrix*2
        for k = 1:3
		    % определяем максимальное значение сумматоров
		    if adaptive_filter_struct_max.Adaptive_filter_sum_array_max(n,k) < adaptive_filter_structure.sum_array_out(n,k)
			    adaptive_filter_struct_max.Adaptive_filter_sum_array_max(n,k) = adaptive_filter_structure.sum_array_out(n,k);
            end
            % определяем максимальную разрядность сумматоров
		    if adaptive_filter_struct_max.Adaptive_filter_sum_total_width(n,k) < adaptive_filter_structure.sum_int_width_total(n,k) 
			    adaptive_filter_struct_max.Adaptive_filter_sum_total_width(n,k) = adaptive_filter_structure.sum_int_width_total(n,k);
            end
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
% title('Относительная ошибка между определителями 2x2, найденных с помощью прямого нахождения в double vs Функции Матлаб')
% ylabel('Величина ошибки') 
% xlabel('Номер отсчета') 
% subplot(2,1,2)
% plot(relativeError_DetM_2x2_Myfunc_double_vs_Myfunc_int, '-o');
% title('Относительная ошибка между определителями 2x2, найденных с помощью прямого нахождения double vs integer')
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
% title('Относительная ошибка между определителями 3x3, найденных с помощью прямого нахождения в double vs Функции Матлаб');
% ylabel('Величина ошибки'); 
% xlabel('Номер отсчета'); 
% subplot(2,1,2)
% plot(relativeError_DetM_3x3_Myfunc_double_vs_Myfunc_int, '-o');
% title('Относительная ошибка между определителями 3x3, найденных с помощью прямого нахождения double vs integer');
% ylabel('Величина ошибки');
% xlabel('Номер отсчета'); 
	
%% 4x4
% relativeError_DetM_4x4_Myfunc_double_vs_Myfunc_int = DetM_4x4_array_dd./double(DetM_4x4_array_int);
% relativeError_DetM_4x4_Myfunc_double_vs_Matlab_LU = Det_4x4_LU_matlab_array./double(DetM_4x4_array_int);
% 
% figure(20)
% subplot(2,1,1)
% plot(relativeError_DetM_4x4_Myfunc_double_vs_Matlab_LU, '-o');
% title('Относительная ошибка между определителями 4x4, найденных с помощью прямого нахождения в double vs Функции Матлаб');
% ylabel('Величина ошибки'); 
% xlabel('Номер отсчета'); 
% subplot(2,1,2)
% plot(relativeError_DetM_4x4_Myfunc_double_vs_Myfunc_int, '-o');
% title('Относительная ошибка между определителями 4x4, найденных с помощью прямого нахождения double vs integer');
% ylabel('Величина ошибки'); 
% xlabel('Номер отсчета'); 

%% 5x5
% relativeError_DetM_5x5_LU_vs_Myfunc = det_x3_shift./double(det_out_shift_int);
% relativeError_DetM_5x5_Myfunc_double_vs_Myfunc_int = det_out_shift./double(det_out_shift_int);
% 
% figure(21)
% subplot(2,1,1)
% plot(relativeError_DetM_5x5_LU_vs_Myfunc, '-o');
% title('Относительная ошибка между определителями, найденных с помощью LU-преобразования и прямого нахождения')
% ylabel('Величина ошибки') 
% xlabel('Номер отсчета') 
% subplot(2,1,2)
% plot(relativeError_DetM_5x5_Myfunc_double_vs_Myfunc_int, '-o');
% title('Относительная ошибка между определителями, найденных с помощью прямого нахождения Double vs Integer')
% ylabel('Величина ошибки') 
% xlabel('Номер отсчета') 
% % 
% relativeError_y_out_matlab_vs_y_out_int = y_array ./ double(y_array_int);
% relativeError_y_out_double_vs_y_out_int = y_array_double ./ double(y_array_int);
% 
% figure(29); 
% subplot(2,1,1)
% plot(relativeError_y_out_matlab_vs_y_out_int);
% title('Относительная ошибка выходного сигнала алгоритма, построенного с помощью матлаб функций и собственных функций в int')
% ylabel('Величина ошибки') 
% xlabel('Номер отсчета')
% subplot(2,1,2)
% plot(relativeError_y_out_double_vs_y_out_int);
% title('Относительная ошибка выходного сигнала алгоритма, построенного с помощью собственных функций в double и int')
% ylabel('Величина ошибки') 
% xlabel('Номер отсчета')

end

function determinate_struct = compare_determinante(s, determinate_struct, sim_options); 
%% определяем макс. значения в определителе 2x2
	for n = 1:sim_options.num_det2x2*2
		% определяем максимальное значение на каждом из 20 умножителей
		if determinate_struct.DetM_2x2_multiplier_total_abs_max(n) < s.Det2x2_mult_abs(n) 
			determinate_struct.DetM_2x2_multiplier_total_abs_max(n) = s.Det2x2_mult_abs(n);
		end
	end 
	for n = 1:sim_options.num_det2x2
		% определяем максимальное значение на каждом из 10 сумматоров
		if determinate_struct.Det2x2_sum_abs_max(n) < s.Det2x2_sum_abs(n) 
			determinate_struct.Det2x2_sum_abs_max(n) = s.Det2x2_sum_abs(n);
		end
	end 
	%% 																		3x3
	for n = 1:sim_options.num_det2x2*3
		% определяем максимальное значение на каждом из 30 умножителей
		if determinate_struct.Mult_DetM_3x3_array_max(n) < s.Mult_DetM_3x3_array(n) 
			determinate_struct.Mult_DetM_3x3_array_max(n) = s.Mult_DetM_3x3_array(n);
		end		
		% определяем макс. разрядность на каждом из 30 умножителей
		if determinate_struct.Mult_DetM_3x3_array_mult_total_width_max(n) < s.Mult_DetM_3x3_array_mult_total_width(n) 
			determinate_struct.Mult_DetM_3x3_array_mult_total_width_max(n) = s.Mult_DetM_3x3_array_mult_total_width(n);
		end
	end
	for n = 1:sim_options.num_det2x2
		% определяем максимальное значение на каждом из 10 пресумматорах
		if determinate_struct.DetM_3x3_int_pre_sum_array_max(n) < s.DetM_3x3_int_pre_sum_array(n) 
			determinate_struct.DetM_3x3_int_pre_sum_array_max(n) = s.DetM_3x3_int_pre_sum_array(n);
		end		
		% определяем макс. разрядность на каждом из 10 пресумматорах
		if determinate_struct.DetM_3x3_int_pre_sum_width_total_max(n) < s.DetM_3x3_int_pre_sum_width_total(n) 
			determinate_struct.DetM_3x3_int_pre_sum_width_total_max(n) = s.DetM_3x3_int_pre_sum_width_total(n);
		end
		% определяем максимальное значение на каждом из 10 сумматорах
		if determinate_struct.DetM_3x3_int_sum_array_max(n) < s.DetM_3x3_int_sum_array(n) 
			determinate_struct.DetM_3x3_int_sum_array_max(n) = s.DetM_3x3_int_sum_array(n);
		end		
		% определяем макс. разрядность на каждом из 10 сумматорах
		if determinate_struct.DetM_3x3_int_sum_width_total_max(n) < s.DetM_3x3_int_sum_width_total(n) 
			determinate_struct.DetM_3x3_int_sum_width_total_max(n) = s.DetM_3x3_int_sum_width_total(n);
		end
	end
	%% 																		4x4 
	for n = 1:sim_options.num_det2x2*2
		% определяем максимальное значение на каждом из 20 умножителей
		if  determinate_struct.DetM_4x4_int_mult_array_max(n) < s.DetM_4x4_int_mult_array(n) 
			determinate_struct.DetM_4x4_int_mult_array_max(n) = s.DetM_4x4_int_mult_array(n);
		end		
		% определяем макс. разрядность на каждом из 20 умножителей
		if  determinate_struct.DetM_4x4_int_mult_width_total_max(n) < s.DetM_4x4_int_mult_width_total(n) 
			determinate_struct.DetM_4x4_int_mult_width_total_max(n) = s.DetM_4x4_int_mult_width_total(n);
		end
	end
	for n = 1:sim_options.num_det2x2
		% определяем максимальное значение на каждом из 10 пресумматоров
		if  determinate_struct.DetM_4x4_int_pre_sum_array_max(n) < s.DetM_4x4_int_pre_sum_array(n) 
			determinate_struct.DetM_4x4_int_pre_sum_array_max(n) = s.DetM_4x4_int_pre_sum_array(n);
		end		
		% определяем макс. разрядность на каждом из 10 пресумматоров
		if  determinate_struct.DetM_4x4_int_pre_sum_width_total_max(n) < s.DetM_4x4_int_pre_sum_width_total(n) 
			determinate_struct.DetM_4x4_int_pre_sum_width_total_max(n) = s.DetM_4x4_int_pre_sum_width_total(n);
		end
	end
	for n = 1:5
		% определяем максимальное значение на каждом из 5 сумматоров
		if  determinate_struct.DetM_4x4_int_sum_array_max(n) < s.DetM_4x4_int_sum_array(n) 
			determinate_struct.DetM_4x4_int_sum_array_max(n) = s.DetM_4x4_int_sum_array(n);
		end		
		% определяем макс. разрядность на каждом из 5 сумматоров
		if  determinate_struct.DetM_4x4_int_sum_array_width_total_max(n) < s.DetM_4x4_int_sum_array_width_total(n) 
			determinate_struct.DetM_4x4_int_sum_array_width_total_max(n) = s.DetM_4x4_int_sum_array_width_total(n);
		end
	end
	%% 																		5x5
	for n = 1:5
		% определяем максимальное значение на каждом из 5 сумматоров
		if  determinate_struct.DetM_5x5_int_mult_array_max(n) < s.DetM_5x5_int_mult_array(n) 
			determinate_struct.DetM_5x5_int_mult_array_max(n) = s.DetM_5x5_int_mult_array(n);
		end		
		% определяем макс. разрядность на каждом из 5 сумматоров
		if  determinate_struct.DetM_5x5_int_mult_array_width_total_max(n) < s.DetM_5x5_int_mult_array_width_total(n) 
			determinate_struct.DetM_5x5_int_mult_array_width_total_max(n) = s.DetM_5x5_int_mult_array_width_total(n);
		end
	end
	for n = 1:2
		% определяем максимальное значение на каждом из 5 сумматоров
		if  determinate_struct.DetM_5x5_int_pre_sum1_array_max(n) < s.DetM_5x5_int_pre_sum1_array(n) 
			determinate_struct.DetM_5x5_int_pre_sum1_array_max(n) = s.DetM_5x5_int_pre_sum1_array(n);
		end		
		% определяем макс. разрядность на каждом из 5 сумматоров
		if  determinate_struct.DetM_5x5_int_pre_sum1_array_width_total_max(n) < s.DetM_5x5_int_pre_sum1_array_width_total(n) 
			determinate_struct.DetM_5x5_int_pre_sum1_array_width_total_max(n) = s.DetM_5x5_int_pre_sum1_array_width_total(n);
		end
	end
	% пресумматор2 определителя 5х5
	if  determinate_struct.DetM_5x5_int_sum3_abs_max < s.DetM_5x5_int_sum3_abs 
		determinate_struct.DetM_5x5_int_sum3_abs_max = s.DetM_5x5_int_sum3_abs;
	end
	% разрядность пресумматора2 определителя 5х5
	if  determinate_struct.DetM_5x5_int_sum3_width_total_max < s.DetM_5x5_int_sum3_width_total 
		determinate_struct.DetM_5x5_int_sum3_width_total_max = s.DetM_5x5_int_sum3_width_total;
	end
	% сумматор определителя 5х5
	if  determinate_struct.DetM_5x5_int_abs_max < s.DetM_5x5_int_abs 
		determinate_struct.DetM_5x5_int_abs_max = s.DetM_5x5_int_abs;
	end
	% разрядность сумматора определителя 5х5
	if  determinate_struct.DetM_5x5_int_width_total_max < s.DetM_5x5_int_width_total 
		determinate_struct.DetM_5x5_int_width_total_max = s.DetM_5x5_int_width_total;
    end
    % начальный определитель
    if determinate_struct.Det_x3_int_max < s.DetM_5x5_int_abs
        determinate_struct.Det_x3_int_max = s.DetM_5x5_int_abs;
    end
end