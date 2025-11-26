
function [DetM_5x5, DetM_5x5_int, DetM_2x2_abs, Det_2x2_LU_matlab, DetM_3x3_array, Det_3x3_LU_matlab, DetM_4x4_array, Det_4x4_LU_matlab, ...
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
] = determinate(data_in, data_in_int, sim_options)

a = data_in;
a_int = data_in_int;

%% Находим матрицы 2x2
s = struct;
e = 0;
t = 4;
for n = 1:4
    for i = 1:t
        for j = 1:1
            s.a{i+e}(:,j) = a(4:end,n);
            s.a_int{i+e}(:,j) = a_int(4:end,n);
        end
        s.a{i+e}(:,2) = a(4:end,n+i);
        s.a_int{i+e}(:,2) = a_int(4:end,n+i);
    end
    e = e+i;
    t = t-1;
end

DetM_2x2 = zeros(sim_options.num_det2x2,1);
DetM_2x2_int = cast(zeros(sim_options.num_det2x2,1), sim_options.type_2x2_det);

Det_2x2_LU_matlab = zeros(sim_options.num_det2x2,1);

Det2x2_mult1_abs = cast(zeros(sim_options.num_det2x2,1), sim_options.type_2x2_det);
Det2x2_mult2_abs = cast(zeros(sim_options.num_det2x2,1), sim_options.type_2x2_det);
Det2x2_sum_abs = cast(zeros(sim_options.num_det2x2,1), 	 sim_options.type_2x2_det);
Det2x2_mult1_abs_max = cast(zeros(sim_options.num_det2x2,1), sim_options.type_2x2_det);
Det2x2_mult2_abs_max = cast(zeros(sim_options.num_det2x2,1), sim_options.type_2x2_det);
Det2x2_sum_abs_max = cast(zeros(sim_options.num_det2x2,1), 	 sim_options.type_2x2_det);
DetM_2x2_multiplier_total_abs = cast(zeros(sim_options.num_det2x2*2,1), sim_options.type_2x2_det);

%% Определитель 3х3
Mult_DetM_3x3_n_1_abs_max = cast(zeros(sim_options.num_det2x2*3,1), sim_options.type_3x3_det);
Mult_DetM_3x3_n_1_total_width_max = cast(zeros(sim_options.num_det2x2*3,1), sim_options.type_3x3_det);
Mult_DetM_3x3_n_2_abs_max = cast(zeros(sim_options.num_det2x2*3,1), sim_options.type_3x3_det);
Mult_DetM_3x3_n_2_total_width_max = cast(zeros(sim_options.num_det2x2*3,1), sim_options.type_3x3_det);
Mult_DetM_3x3_n_3_abs_max = cast(zeros(sim_options.num_det2x2*3,1), sim_options.type_3x3_det);
Mult_DetM_3x3_n_3_total_width_max = cast(zeros(sim_options.num_det2x2*3,1), sim_options.type_3x3_det);
Mult_DetM_3x3_n_4_abs_max = cast(zeros(sim_options.num_det2x2*3,1), sim_options.type_3x3_det);
Mult_DetM_3x3_n_4_total_width_max = cast(zeros(sim_options.num_det2x2*3,1), sim_options.type_3x3_det);
Mult_DetM_3x3_n_5_abs_max = cast(zeros(sim_options.num_det2x2*3,1), sim_options.type_3x3_det);
Mult_DetM_3x3_n_5_total_width_max = cast(zeros(sim_options.num_det2x2*3,1), sim_options.type_3x3_det);
Mult_DetM_3x3_n_6_abs_max = cast(zeros(sim_options.num_det2x2*3,1), sim_options.type_3x3_det);
Mult_DetM_3x3_n_6_total_width_max = cast(zeros(sim_options.num_det2x2*3,1), sim_options.type_3x3_det);
Mult_DetM_3x3_n_7_abs_max = cast(zeros(sim_options.num_det2x2*3,1), sim_options.type_3x3_det);
Mult_DetM_3x3_n_7_total_width_max = cast(zeros(sim_options.num_det2x2*3,1), sim_options.type_3x3_det);
Mult_DetM_3x3_n_8_abs_max = cast(zeros(sim_options.num_det2x2*3,1), sim_options.type_3x3_det);
Mult_DetM_3x3_n_8_total_width_max = cast(zeros(sim_options.num_det2x2*3,1), sim_options.type_3x3_det);
Mult_DetM_3x3_n_9_abs_max = cast(zeros(sim_options.num_det2x2*3,1), sim_options.type_3x3_det);
Mult_DetM_3x3_n_9_total_width_max = cast(zeros(sim_options.num_det2x2*3,1), sim_options.type_3x3_det);
Mult_DetM_3x3_n_10_abs_max = cast(zeros(sim_options.num_det2x2*3,1), sim_options.type_3x3_det);
Mult_DetM_3x3_n_10_total_width_max = cast(zeros(sim_options.num_det2x2*3,1), sim_options.type_3x3_det);
DetM_3x3_n_44_int_sum1_abs_max = cast(zeros(1,1), sim_options.type_3x3_det);
DetM_3x3_n_44_int_sum1_width_total_max = cast(zeros(1,1), sim_options.type_3x3_det);
DetM_3x3_n_34_int_sum1_abs_max = cast(zeros(1,1), sim_options.type_3x3_det);
DetM_3x3_n_34_int_sum1_width_total_max = cast(zeros(1,1), sim_options.type_3x3_det);
DetM_3x3_n_33_int_sum1_abs_max = cast(zeros(1,1), sim_options.type_3x3_det);
DetM_3x3_n_33_int_sum1_width_total_max = cast(zeros(1,1), sim_options.type_3x3_det);
DetM_3x3_n_24_int_sum1_abs_max = cast(zeros(1,1), sim_options.type_3x3_det);
DetM_3x3_n_24_int_sum1_width_total_max = cast(zeros(1,1), sim_options.type_3x3_det);
DetM_3x3_n_23_int_sum1_abs_max = cast(zeros(1,1), sim_options.type_3x3_det);
DetM_3x3_n_23_int_sum1_width_total_max = cast(zeros(1,1), sim_options.type_3x3_det);
DetM_3x3_n_22_int_sum1_abs_max = cast(zeros(1,1), sim_options.type_3x3_det);
DetM_3x3_n_22_int_sum1_width_total_max = cast(zeros(1,1), sim_options.type_3x3_det);
DetM_3x3_n_14_int_sum1_abs_max = cast(zeros(1,1), sim_options.type_3x3_det);
DetM_3x3_n_14_int_sum1_width_total_max = cast(zeros(1,1), sim_options.type_3x3_det);
DetM_3x3_n_13_int_sum1_abs_max = cast(zeros(1,1), sim_options.type_3x3_det);
DetM_3x3_n_13_int_sum1_width_total_max = cast(zeros(1,1), sim_options.type_3x3_det);
DetM_3x3_n_12_int_sum1_abs_max = cast(zeros(1,1), sim_options.type_3x3_det);
DetM_3x3_n_12_int_sum1_width_total_max = cast(zeros(1,1), sim_options.type_3x3_det);
DetM_3x3_n_11_int_sum1_abs_max = cast(zeros(1,1), sim_options.type_3x3_det);
DetM_3x3_n_11_int_sum1_width_total_max = cast(zeros(1,1), sim_options.type_3x3_det);


width_total_mult1_max = cast(zeros(sim_options.num_det2x2,1), 	 sim_options.type_2x2_det);
width_total_mult2_max = cast(zeros(sim_options.num_det2x2,1), 	 sim_options.type_2x2_det);
width_total_sum_max = cast(zeros(sim_options.num_det2x2,1), 	 sim_options.type_2x2_det);
%%
Det_3x3_LU_matlab 						= 	zeros(sim_options.num_det2x2,1);
Mult_DetM_3x3_array 					= 	cast(zeros(sim_options.num_det2x2*3,1), sim_options.type_3x3_det);
Mult_DetM_3x3_array_mult_total_width 	= 	cast(zeros(sim_options.num_det2x2*3,1), sim_options.type_3x3_det);
DetM_3x3_int_pre_sum_array 				= 	cast(zeros(sim_options.num_det2x2,1), sim_options.type_3x3_det);
DetM_3x3_int_pre_sum_width_total 		= 	cast(zeros(sim_options.num_det2x2,1), sim_options.type_3x3_det);
DetM_3x3_int_sum_array 					= 	cast(zeros(sim_options.num_det2x2,1), sim_options.type_3x3_det);
DetM_3x3_int_sum_width_total 			=	cast(zeros(sim_options.num_det2x2,1), sim_options.type_3x3_det);

DetM_4x4_int_mult_array	 				= cast(zeros(sim_options.num_det2x2*2,1), sim_options.type_4x4_det);
DetM_4x4_int_mult_width_total 			= cast(zeros(sim_options.num_det2x2*2,1), sim_options.type_4x4_det);
DetM_4x4_int_pre_sum_array 				= cast(zeros(sim_options.num_det2x2,1), sim_options.type_4x4_det);
DetM_4x4_int_pre_sum_width_total 		= cast(zeros(sim_options.num_det2x2,1), sim_options.type_4x4_det);
DetM_4x4_int_sum_array 					= cast(zeros(5,1), sim_options.type_4x4_det);
DetM_4x4_int_sum_array_width_total 		= cast(zeros(5,1), sim_options.type_4x4_det);

DetM_5x5_int_mult_array 				= cast(zeros(5,1), sim_options.type_5x5_det);
DetM_5x5_int_mult_array_width_total 	= cast(zeros(5,1), sim_options.type_5x5_det);
DetM_5x5_int_pre_sum1_array 			= cast(zeros(2,1), sim_options.type_5x5_det);
DetM_5x5_int_pre_sum1_array_width_total = cast(zeros(2,1), sim_options.type_5x5_det);

% if sim_options.enable_mask == true
%     width_mult1_det_2x2 = readmatrix(width_mult_txt);
%     width_mult2_det_2x2 = readmatrix(width_sum_txt);
%     width_sum_det_2x2 = readmatrix(width_sum_txt);
% end

%% Находим определители матриц 2x2
for i = 1:sim_options.num_det2x2

    ee = det(s.a{i});
    if (ee < 0)
        Det_2x2_LU_matlab(i) = ee * -1;
    else
        Det_2x2_LU_matlab(i) = ee;
    end
    %% Определитель 2х2

    DetM_2x2(i) = s.a{i}(1,1) * s.a{i}(2,2) - s.a{i}(2,1) * s.a{i}(1,2);
	
    %% 
    [mult1_int(i), mult1_overflow(i), Det2x2_mult1_abs(i), width_total_mult1(i)] = mult(cast(s.a_int{i}(1,1),sim_options.type_2x2_det), ...
        cast(s.a_int{i}(2,2),sim_options.type_2x2_det), sim_options.type_2x2_det, sim_options.width_fractional);

     % Наложение маски на первый умножитель
     if sim_options.enable_mask == true
        c1 = bitmask(mult1_int(i), sim_options.type_2x2_det, width_mult1(i));
        if c1 ~= mult1_int(i)
            disp('Bit mask error mult 1 det 2x2');
            % c2 = bitmask(sum(i+1,n), sim_options.int_size, width_sum(i+1));
            disp(width_total_mult1(i));
            disp({c1, mult1_int(i)});
            disp({sim_options.freq, sim_options.SNR});
        end
        mult1_int(i) = c1;
    else
        % записываем макс значение на сумматоре
        if Det2x2_mult1_abs_max(i) < Det2x2_mult1_abs(i)
            Det2x2_mult1_abs_max(i) = Det2x2_mult1_abs(i);
        end

        % записываем макс значение суммы разрядностей сумматоров
        if width_total_mult1_max(i) < width_total_mult1(i)
            width_total_mult1_max(i) = width_total_mult1(i);
        end
    end

    %%
    [mult2_int(i), mult2_overflow(i), Det2x2_mult2_abs(i), width_total_mult2(i)] = mult(cast(s.a_int{i}(2,1),sim_options.type_2x2_det), ...
        cast(s.a_int{i}(1,2),sim_options.type_2x2_det), sim_options.type_2x2_det, sim_options.width_fractional);

     % Наложение маски на второй умножитель
     if sim_options.enable_mask == true
        c1 = bitmask(mult2_int(i), sim_options.type_2x2_det, width_mult2(i));
        if c1 ~= mult2_int(i)
            disp('Bit mask error mult 2 det 2x2');
            disp(width_total_mult2(i));
            disp({c1, mult2_int(i)});
            disp({sim_options.freq, sim_options.SNR});
        end
        mult2_int(i) = c1;
    else
        % записываем макс значение на сумматоре
        if Det2x2_mult2_abs_max(i) < Det2x2_mult2_abs(i)
            Det2x2_mult2_abs_max(i) = Det2x2_mult2_abs(i);
        end

        % записываем макс значение суммы разрядностей сумматоров
        if width_total_mult2_max(i) < width_total_mult2(i)
            width_total_mult2_max(i) = width_total_mult2(i);
        end
     end

	[DetM_2x2_int(i), sum_overflow(i), Det2x2_sum_abs(i), width_total_sum(i)] = adder(cast(mult1_int(i),sim_options.type_2x2_det), ...
        cast(-mult2_int(i),sim_options.type_2x2_det), sim_options.type_2x2_det, sim_options.width_fractional);

     % Наложение маски на сумматор
     if sim_options.enable_mask == true
        c1 = bitmask(DetM_2x2_int(i), sim_options.type_2x2_det, width_sum(i));
        if c1 ~= DetM_2x2_int(i)
            disp('Bit mask error sum det 2x2');
            % c2 = bitmask(sum(i+1,n), sim_options.int_size, width_sum(i+1));
            disp(width_total_sum(i));
            disp({c1, DetM_2x2_int(i)});
            disp({sim_options.freq, sim_options.SNR});
        end
        DetM_2x2_int(i) = c1;
    else
        % записываем макс значение на сумматоре
        if Det2x2_sum_abs_max(i) < Det2x2_sum_abs(i)
            Det2x2_sum_abs_max(i) = Det2x2_sum_abs(i);
        end

        % записываем макс значение суммы разрядностей сумматоров
        if width_total_sum_max(i) < width_total_sum(i)
            width_total_sum_max(i) = width_total_sum(i);
        end
     end

    %% Проверка переполнения умножителя
    if (mult1_overflow(i) == 1 || mult2_overflow(i) == 1)
		disp('Mult in function det2x2 overflow');
        disp(sim_options.width_hilbert);
        disp({s.a_int{i}});
    end
    % 
	%% Проверка выходной разрядности умножителя
    if (width_total_mult1(i) > sim_options.width_hilbert || width_total_mult2(i) > sim_options.width_hilbert)
		disp('Mult total width in function det2x2 higher than 64');
        disp(sim_options.width_hilbert);
        disp({s.a_int{i}});
    end
    % 
	%% Проверка переполнения сумматора
    if (sum_overflow(i) == 1)
		disp('Adder in function det2x2 overflow');
        disp(sim_options.width_hilbert);
        disp({s.a_int{i}});
    end
    % 
	%% Проверка выходной разрядности сумматора
    if (width_total_sum(i) > sim_options.width_hilbert)
		disp('Sum total width in function det2x2 higher than 64');
        disp(sim_options.width_hilbert);
        disp({width_total_sum(i)});
    end
end


%% 3x3
% double
Mult_DetM_2x2_n_7_a31 = a(3,1) * DetM_2x2(7);
Mult_DetM_2x2_n_7_a33 = a(3,3) * DetM_2x2(7);
Mult_DetM_2x2_n_7_a34 = a(3,4) * DetM_2x2(7);

Mult_DetM_2x2_n_9_a31 = a(3,1) * DetM_2x2(9);
Mult_DetM_2x2_n_9_a32 = a(3,2) * DetM_2x2(9);
Mult_DetM_2x2_n_9_a34 = a(3,4) * DetM_2x2(9);

Mult_DetM_2x2_n_8_a31 = a(3,1) * DetM_2x2(8);
Mult_DetM_2x2_n_8_a32 = a(3,2) * DetM_2x2(8);
Mult_DetM_2x2_n_8_a35 = a(3,5) * DetM_2x2(8);

Mult_DetM_2x2_n_6_a31 = a(3,1) * DetM_2x2(6);
Mult_DetM_2x2_n_6_a33 = a(3,3) * DetM_2x2(6);
Mult_DetM_2x2_n_6_a35 = a(3,5) * DetM_2x2(6);

Mult_DetM_2x2_n_5_a31 = a(3,1) * DetM_2x2(5);
Mult_DetM_2x2_n_5_a34 = a(3,4) * DetM_2x2(5);
Mult_DetM_2x2_n_5_a35 = a(3,5) * DetM_2x2(5);

Mult_DetM_2x2_n_4_a32 = a(3,2) * DetM_2x2(4);
Mult_DetM_2x2_n_4_a33 = a(3,3) * DetM_2x2(4);
Mult_DetM_2x2_n_4_a34 = a(3,4) * DetM_2x2(4);

Mult_DetM_2x2_n_3_a32 = a(3,2) * DetM_2x2(3);
Mult_DetM_2x2_n_3_a33 = a(3,3) * DetM_2x2(3);
Mult_DetM_2x2_n_3_a35 = a(3,5) * DetM_2x2(3);

Mult_DetM_2x2_n_2_a32 = a(3,2) * DetM_2x2(2);
Mult_DetM_2x2_n_2_a34 = a(3,4) * DetM_2x2(2);
Mult_DetM_2x2_n_2_a35 = a(3,5) * DetM_2x2(2);

Mult_DetM_2x2_n_1_a33 = a(3,3) * DetM_2x2(1);
Mult_DetM_2x2_n_1_a34 = a(3,4) * DetM_2x2(1);
Mult_DetM_2x2_n_1_a35 = a(3,5) * DetM_2x2(1);

Mult_DetM_2x2_n_10_a31 = a(3,1) * DetM_2x2(10);
Mult_DetM_2x2_n_10_a32 = a(3,2) * DetM_2x2(10);
Mult_DetM_2x2_n_10_a33 = a(3,3) * DetM_2x2(10);

%% integer      
index1 = [1, 3, 4];
index2 = [1, 2, 4];
index3 = [1, 2, 5];
index4 = [1, 3, 5];
index5 = [1, 4, 5];
index6 = [2, 3, 4];
index7 = [2, 3, 5];
index8 = [2, 4, 5];
index9 = [3, 4, 5];

for i = 1:3
    [Mult_DetM_3x3_n_1_int(i), Mult_DetM_3x3_n_1_overflow(i), Mult_DetM_3x3_n_1_abs(i), Mult_DetM_3x3_n_1_total_width(i)] = mult(cast(a_int(3,index9(i)),sim_options.type_3x3_det), cast(DetM_2x2_int(1),sim_options.type_3x3_det), sim_options.type_3x3_det, sim_options.width_hilbert);
	
     % Наложение маски на первый умножитель
     if sim_options.enable_mask == true
        c1 = bitmask(Mult_DetM_3x3_n_1_int(i), sim_options.type_3x3_det, width_mult1(i));
        if c1 ~= Mult_DetM_3x3_n_1_int(i)
            disp('Bit mask error mult 1 det 2x2');
            % c2 = bitmask(sum(i+1,n), sim_options.int_size, width_sum(i+1));
            disp(width_total_mult1(i));
            disp({c1, Mult_DetM_3x3_n_1_int(i)});
            disp({sim_options.freq, sim_options.SNR});
        end
        Mult_DetM_3x3_n_1_int(i) = c1;
    else
        % записываем макс значение на сумматоре
        if Mult_DetM_3x3_n_1_abs_max(i) < Mult_DetM_3x3_n_1_abs(i)
            Mult_DetM_3x3_n_1_abs_max(i) = Mult_DetM_3x3_n_1_abs(i);
        end

        % записываем макс значение суммы разрядностей сумматоров
        if Mult_DetM_3x3_n_1_total_width_max(i) < Mult_DetM_3x3_n_1_total_width(i)
            Mult_DetM_3x3_n_1_total_width_max(i) = Mult_DetM_3x3_n_1_total_width(i);
        end
     end
    %% 3x3
    [Mult_DetM_3x3_n_2_int(i), Mult_DetM_3x3_n_2_overflow(i), Mult_DetM_3x3_n_2_abs(i), Mult_DetM_3x3_n_2_total_width(i)] = mult(cast(a_int(3,index8(i)),sim_options.type_3x3_det), cast(DetM_2x2_int(2),sim_options.type_3x3_det), sim_options.type_3x3_det, sim_options.width_hilbert);
	
     % Наложение маски на первый умножитель
     if sim_options.enable_mask == true
        c1 = bitmask(Mult_DetM_3x3_n_2_int(i), sim_options.type_3x3_det, width_mult1(i));
        if c1 ~= Mult_DetM_3x3_n_2_int(i)
            disp('Bit mask error mult 1 det 2x2');
            % c2 = bitmask(sum(i+1,n), sim_options.int_size, width_sum(i+1));
            disp(width_total_mult1(i));
            disp({c1, Mult_DetM_3x3_n_2_int(i)});
            disp({sim_options.freq, sim_options.SNR});
        end
        Mult_DetM_3x3_n_2_int(i) = c1;
    else
        % записываем макс значение на сумматоре
        if Mult_DetM_3x3_n_3_abs_max(i) < Mult_DetM_3x3_n_2_abs(i)
            Mult_DetM_3x3_n_3_abs_max(i) = Mult_DetM_3x3_n_2_abs(i);
        end

        % записываем макс значение суммы разрядностей сумматоров
        if Mult_DetM_3x3_n_2_total_width_max(i) < Mult_DetM_3x3_n_2_total_width(i)
            Mult_DetM_3x3_n_2_total_width_max(i) = Mult_DetM_3x3_n_2_total_width(i);
        end
     end
    %%
    [Mult_DetM_3x3_n_3_int(i), Mult_DetM_3x3_n_3_overflow(i), Mult_DetM_3x3_n_3_abs(i), Mult_DetM_3x3_n_3_total_width(i)] = mult(cast(a_int(3,index7(i)),sim_options.type_3x3_det), cast(DetM_2x2_int(3),sim_options.type_3x3_det), sim_options.type_3x3_det, sim_options.width_hilbert);
	
    % Наложение маски на первый умножитель
    if sim_options.enable_mask == true
       c1 = bitmask(Mult_DetM_3x3_n_3_int(i), sim_options.type_3x3_det, width_mult1(i));
       if c1 ~= Mult_DetM_3x3_n_3_int(i)
           disp('Bit mask error mult 1 det 2x2');
           % c2 = bitmask(sum(i+1,n), sim_options.int_size, width_sum(i+1));
           disp(width_total_mult1(i));
           disp({c1, Mult_DetM_3x3_n_3_int(i)});
           disp({sim_options.freq, sim_options.SNR});
       end
       Mult_DetM_3x3_n_3_int(i) = c1;
    else
       % записываем макс значение на сумматоре
       if Mult_DetM_3x3_n_3_abs_max(i) < Mult_DetM_3x3_n_3_abs(i)
           Mult_DetM_3x3_n_3_abs_max(i) = Mult_DetM_3x3_n_3_abs(i);
       end
	
       % записываем макс значение суммы разрядностей сумматоров
       if Mult_DetM_3x3_n_3_total_width_max(i) < Mult_DetM_3x3_n_3_total_width(i)
           Mult_DetM_3x3_n_3_total_width_max(i) = Mult_DetM_3x3_n_3_total_width(i);
       end
    end
    %%
    [Mult_DetM_3x3_n_4_int(i), Mult_DetM_3x3_n_4_overflow(i), Mult_DetM_3x3_n_4_abs(i), Mult_DetM_3x3_n_4_total_width(i)] = mult(cast(a_int(3,index6(i)),sim_options.type_3x3_det), cast(DetM_2x2_int(4),sim_options.type_3x3_det), sim_options.type_3x3_det, sim_options.width_hilbert);
	
    % Наложение маски на первый умножитель
    if sim_options.enable_mask == true
       c1 = bitmask(Mult_DetM_3x3_n_4_int(i), sim_options.type_3x3_det, width_mult1(i));
       if c1 ~= Mult_DetM_3x3_n_4_int(i)
           disp('Bit mask error mult 1 det 2x2');
           % c2 = bitmask(sum(i+1,n), sim_options.int_size, width_sum(i+1));
           disp(width_total_mult1(i));
           disp({c1, Mult_DetM_3x3_n_4_int(i)});
           disp({sim_options.freq, sim_options.SNR});
       end
       Mult_DetM_3x3_n_4_int(i) = c1;
    else
       % записываем макс значение на сумматоре
       if Mult_DetM_3x3_n_4_abs_max(i) < Mult_DetM_3x3_n_4_abs(i)
           Mult_DetM_3x3_n_4_abs_max(i) = Mult_DetM_3x3_n_4_abs(i);
       end
	
       % записываем макс значение суммы разрядностей сумматоров
       if Mult_DetM_3x3_n_4_total_width_max(i) < Mult_DetM_3x3_n_4_total_width(i)
           Mult_DetM_3x3_n_4_total_width_max(i) = Mult_DetM_3x3_n_4_total_width(i);
       end
    end
    
	%%
    [Mult_DetM_3x3_n_5_int(i), Mult_DetM_3x3_n_5_overflow(i), Mult_DetM_3x3_n_5_abs(i), Mult_DetM_3x3_n_5_total_width(i)] = mult(cast(a_int(3,index5(i)),sim_options.type_3x3_det), cast(DetM_2x2_int(5),sim_options.type_3x3_det), sim_options.type_3x3_det, sim_options.width_hilbert);
	
	% Наложение маски на первый умножитель
    if sim_options.enable_mask == true
       c1 = bitmask(Mult_DetM_3x3_n_5_int(i), sim_options.type_3x3_det, width_mult1(i));
       if c1 ~= Mult_DetM_3x3_n_5_int(i)
           disp('Bit mask error mult 1 det 2x2');
           % c2 = bitmask(sum(i+1,n), sim_options.int_size, width_sum(i+1));
           disp(width_total_mult1(i));
           disp({c1, Mult_DetM_3x3_n_5_int(i)});
           disp({sim_options.freq, sim_options.SNR});
       end
       Mult_DetM_3x3_n_5_int(i) = c1;
    else
       % записываем макс значение на сумматоре
       if Mult_DetM_3x3_n_5_abs_max(i) < Mult_DetM_3x3_n_5_abs(i)
           Mult_DetM_3x3_n_5_abs_max(i) = Mult_DetM_3x3_n_5_abs(i);
       end
	
       % записываем макс значение суммы разрядностей сумматоров
       if Mult_DetM_3x3_n_5_total_width_max(i) < Mult_DetM_3x3_n_5_total_width(i)
           Mult_DetM_3x3_n_5_total_width_max(i) = Mult_DetM_3x3_n_5_total_width(i);
       end
    end
	
	%%
	[Mult_DetM_3x3_n_6_int(i), Mult_DetM_3x3_n_6_overflow(i), Mult_DetM_3x3_n_6_abs(i), Mult_DetM_3x3_n_6_total_width(i)] = mult(cast(a_int(3,index4(i)),sim_options.type_3x3_det), cast(DetM_2x2_int(6),sim_options.type_3x3_det), sim_options.type_3x3_det, sim_options.width_hilbert);
	
		% Наложение маски на первый умножитель
    if sim_options.enable_mask == true
       c1 = bitmask(Mult_DetM_3x3_n_6_int(i), sim_options.type_3x3_det, width_mult1(i));
       if c1 ~= Mult_DetM_3x3_n_6_int(i)
           disp('Bit mask error mult 1 det 2x2');
           % c2 = bitmask(sum(i+1,n), sim_options.int_size, width_sum(i+1));
           disp(width_total_mult1(i));
           disp({c1, Mult_DetM_3x3_n_6_int(i)});
           disp({sim_options.freq, sim_options.SNR});
       end
       Mult_DetM_3x3_n_6_int(i) = c1;
    else
       % записываем макс значение на сумматоре
       if Mult_DetM_3x3_n_6_abs_max(i) < Mult_DetM_3x3_n_6_abs(i)
           Mult_DetM_3x3_n_6_abs_max(i) = Mult_DetM_3x3_n_6_abs(i);
       end
	
       % записываем макс значение суммы разрядностей сумматоров
       if Mult_DetM_3x3_n_6_total_width_max(i) < Mult_DetM_3x3_n_6_total_width(i)
           Mult_DetM_3x3_n_6_total_width_max(i) = Mult_DetM_3x3_n_6_total_width(i);
       end
    end
	
	%%
	
	[Mult_DetM_3x3_n_7_int(i), Mult_DetM_3x3_n_7_overflow(i), Mult_DetM_3x3_n_7_abs(i), Mult_DetM_3x3_n_7_total_width(i)] = mult(cast(a_int(3,index1(i)),sim_options.type_3x3_det), cast(DetM_2x2_int(7),sim_options.type_3x3_det), sim_options.type_3x3_det, sim_options.width_hilbert);
	
	% Наложение маски на первый умножитель
    if sim_options.enable_mask == true
       c1 = bitmask(Mult_DetM_3x3_n_7_int(i), sim_options.type_3x3_det, width_mult1(i));
       if c1 ~= Mult_DetM_3x3_n_7_int(i)
           disp('Bit mask error mult 1 det 2x2');
           % c2 = bitmask(sum(i+1,n), sim_options.int_size, width_sum(i+1));
           disp(width_total_mult1(i));
           disp({c1, Mult_DetM_3x3_n_7_int(i)});
           disp({sim_options.freq, sim_options.SNR});
       end
       Mult_DetM_3x3_n_7_int(i) = c1;
    else
       % записываем макс значение на сумматоре
       if Mult_DetM_3x3_n_7_abs_max(i) < Mult_DetM_3x3_n_7_abs(i)
           Mult_DetM_3x3_n_7_abs_max(i) = Mult_DetM_3x3_n_7_abs(i);
       end
	
       % записываем макс значение суммы разрядностей сумматоров
       if Mult_DetM_3x3_n_7_total_width_max(i) < Mult_DetM_3x3_n_7_total_width(i)
           Mult_DetM_3x3_n_7_total_width_max(i) = Mult_DetM_3x3_n_7_total_width(i);
       end
    end
	
	%%
	[Mult_DetM_3x3_n_8_int(i), Mult_DetM_3x3_n_8_overflow(i), Mult_DetM_3x3_n_8_abs(i), Mult_DetM_3x3_n_8_total_width(i)] = mult(cast(a_int(3,index3(i)),sim_options.type_3x3_det), cast(DetM_2x2_int(8),sim_options.type_3x3_det), sim_options.type_3x3_det, sim_options.width_hilbert);
	
	% Наложение маски на первый умножитель
    if sim_options.enable_mask == true
       c1 = bitmask(Mult_DetM_3x3_n_8_int(i), sim_options.type_3x3_det, width_mult1(i));
       if c1 ~= Mult_DetM_3x3_n_8_int(i)
           disp('Bit mask error mult 1 det 2x2');
           % c2 = bitmask(sum(i+1,n), sim_options.int_size, width_sum(i+1));
           disp(width_total_mult1(i));
           disp({c1, Mult_DetM_3x3_n_8_int(i)});
           disp({sim_options.freq, sim_options.SNR});
       end
       Mult_DetM_3x3_n_8_int(i) = c1;
    else
       % записываем макс значение на сумматоре
       if Mult_DetM_3x3_n_8_abs_max(i) < Mult_DetM_3x3_n_8_abs(i)
           Mult_DetM_3x3_n_8_abs_max(i) = Mult_DetM_3x3_n_8_abs(i);
       end
	
       % записываем макс значение суммы разрядностей сумматоров
       if Mult_DetM_3x3_n_8_total_width_max(i) < Mult_DetM_3x3_n_8_total_width(i)
           Mult_DetM_3x3_n_8_total_width_max(i) = Mult_DetM_3x3_n_8_total_width(i);
       end
    end
	
	%%
	[Mult_DetM_3x3_n_9_int(i), Mult_DetM_3x3_n_9_overflow(i), Mult_DetM_3x3_n_9_abs(i), Mult_DetM_3x3_n_9_total_width(i)] = mult(cast(a_int(3,index2(i)),sim_options.type_3x3_det), cast(DetM_2x2_int(9),sim_options.type_3x3_det), sim_options.type_3x3_det, sim_options.width_hilbert);
	
	% Наложение маски на первый умножитель
    if sim_options.enable_mask == true
       c1 = bitmask(Mult_DetM_3x3_n_9_int(i), sim_options.type_3x3_det, width_mult1(i));
       if c1 ~= Mult_DetM_3x3_n_9_int(i)
           disp('Bit mask error mult 1 det 2x2');
           % c2 = bitmask(sum(i+1,n), sim_options.int_size, width_sum(i+1));
           disp(width_total_mult1(i));
           disp({c1, Mult_DetM_3x3_n_9_int(i)});
           disp({sim_options.freq, sim_options.SNR});
       end
       Mult_DetM_3x3_n_9_int(i) = c1;
    else
       % записываем макс значение на сумматоре
       if Mult_DetM_3x3_n_9_abs_max(i) < Mult_DetM_3x3_n_9_abs(i)
           Mult_DetM_3x3_n_9_abs_max(i) = Mult_DetM_3x3_n_9_abs(i);
       end
	
       % записываем макс значение суммы разрядностей сумматоров
       if Mult_DetM_3x3_n_9_total_width_max(i) < Mult_DetM_3x3_n_9_total_width(i)
           Mult_DetM_3x3_n_9_total_width_max(i) = Mult_DetM_3x3_n_9_total_width(i);
       end
    end
	
	%%
	[Mult_DetM_3x3_n_10_int(i), Mult_DetM_3x3_n_10_overflow(i), Mult_DetM_3x3_n_10_abs(i), Mult_DetM_3x3_n_10_total_width(i)] = mult(cast(a_int(3,i),sim_options.type_3x3_det), 	cast(DetM_2x2_int(10),sim_options.type_3x3_det),sim_options.type_3x3_det, sim_options.width_hilbert);
	
	% Наложение маски на первый умножитель
    if sim_options.enable_mask == true
       c1 = bitmask(Mult_DetM_3x3_n_10_int(i), sim_options.type_3x3_det, width_mult1(i));
       if c1 ~= Mult_DetM_3x3_n_10_int(i)
           disp('Bit mask error mult 1 det 2x2');
           % c2 = bitmask(sum(i+1,n), sim_options.int_size, width_sum(i+1));
           disp(width_total_mult1(i));
           disp({c1, Mult_DetM_3x3_n_10_int(i)});
           disp({sim_options.freq, sim_options.SNR});
       end
       Mult_DetM_3x3_n_10_int(i) = c1;
    else
       % записываем макс значение на сумматоре
       if Mult_DetM_3x3_n_10_abs_max(i) < Mult_DetM_3x3_n_10_abs(i)
           Mult_DetM_3x3_n_10_abs_max(i) = Mult_DetM_3x3_n_10_abs(i);
       end
	
       % записываем макс значение суммы разрядностей сумматоров
       if Mult_DetM_3x3_n_10_total_width_max(i) < Mult_DetM_3x3_n_10_total_width(i)
           Mult_DetM_3x3_n_10_total_width_max(i) = Mult_DetM_3x3_n_10_total_width(i);
       end
    end
	
    %% Проверка переполнения умножителя
    if (Mult_DetM_3x3_n_1_overflow(i) == 1)
		disp('Mult in function det3x3 1 overflow');
        disp({Mult_DetM_3x3_n_1_int(i)});
    end
	%% Проверка выходной разрядности умножителя
    if (Mult_DetM_3x3_n_1_total_width(i) > sim_options.width_hilbert-1)
		disp('Mult total width in function det3x3 1 higher than 64');
        disp({Mult_DetM_3x3_n_1_int(i)});
    end
    %% Проверка переполнения умножителя
    if (Mult_DetM_3x3_n_2_overflow(i) == 1)
		disp('Mult in function det3x3 2 overflow');
        disp({Mult_DetM_3x3_n_2_int(i)});
    end
	%% Проверка выходной разрядности умножителя
    if (Mult_DetM_3x3_n_2_total_width(i) > sim_options.width_hilbert-1)
		disp('Mult total width in function det3x3 2 higher than 64');
        disp({Mult_DetM_3x3_n_2_int(i)});
    end
    %% Проверка переполнения умножителя
    if (Mult_DetM_3x3_n_3_overflow(i) == 1)
		disp('Mult in function det3x3 3 overflow');
        disp({Mult_DetM_3x3_n_3_int(i)});
    end
	%% Проверка выходной разрядности умножителя
    if (Mult_DetM_3x3_n_3_total_width(i) > sim_options.width_hilbert-1)
		disp('Mult total width in function det3x3 3 higher than 64');
        disp({Mult_DetM_3x3_n_3_int(i)});
    end
    %% Проверка переполнения умножителя
    if (Mult_DetM_3x3_n_4_overflow(i) == 1)
		disp('Mult in function det3x3 4 overflow');
        disp({Mult_DetM_3x3_n_4_int(i)});
    end
	%% Проверка выходной разрядности умножителя
    if (Mult_DetM_3x3_n_4_total_width(i) > sim_options.width_hilbert-1)
		disp('Mult total width in function det3x3 4 higher than 64');
        disp({Mult_DetM_3x3_n_4_int(i)});
    end
    %% Проверка переполнения умножителя
    if (Mult_DetM_3x3_n_5_overflow(i) == 1)
		disp('Mult in function det3x3 5 overflow');
        disp({Mult_DetM_3x3_n_5_int(i)});
    end
	%% Проверка выходной разрядности умножителя
    if (Mult_DetM_3x3_n_5_total_width(i) > sim_options.width_hilbert-1)
		disp('Mult total width in function det3x3 5 higher than 64');
        disp({Mult_DetM_3x3_n_5_int(i)});
    end
    %% Проверка переполнения умножителя
    if (Mult_DetM_3x3_n_6_overflow(i) == 1)
		disp('Mult in function det3x3 6 overflow');
        disp({Mult_DetM_3x3_n_6_int(i)});
    end
	%% Проверка выходной разрядности умножителя
    if (Mult_DetM_3x3_n_6_total_width(i) > sim_options.width_hilbert-1)
		disp('Mult total width in function det3x3 6 higher than 64');
        disp({Mult_DetM_3x3_n_6_int(i)});
    end
    %% Проверка переполнения умножителя
    if (Mult_DetM_3x3_n_7_overflow(i) == 1)
		disp('Mult in function det3x3 7 overflow');
        disp({Mult_DetM_3x3_n_7_int(i)});
    end
	%% Проверка выходной разрядности умножителя
    if (Mult_DetM_3x3_n_7_total_width(i) > sim_options.width_hilbert-1)
		disp('Mult total width in function det3x3 7 higher than 64');
        disp({Mult_DetM_3x3_n_7_int(i)});
    end
    %% Проверка переполнения умножителя
    if (Mult_DetM_3x3_n_8_overflow(i) == 1)
		disp('Mult in function det3x3 8 overflow');
        disp({Mult_DetM_3x3_n_8_int(i)});
    end
	%% Проверка выходной разрядности умножителя
    if (Mult_DetM_3x3_n_8_total_width(i) > sim_options.width_hilbert-1)
		disp('Mult total width in function det3x3 8 higher than 64');
        disp({Mult_DetM_3x3_n_8_int(i)});
    end
    %% Проверка переполнения умножителя
    if (Mult_DetM_3x3_n_9_overflow(i) == 1)
		disp('Mult in function det3x3 9 overflow');
        disp({Mult_DetM_3x3_n_9_int(i)});
    end
	%% Проверка выходной разрядности умножителя
    if (Mult_DetM_3x3_n_9_total_width(i) > sim_options.width_hilbert-1)
		disp('Mult total width in function det3x3 9 higher than 64');
        disp({Mult_DetM_3x3_n_9_int(i)});
    end
    %% Проверка переполнения умножителя
    if (Mult_DetM_3x3_n_10_overflow(i) == 1)
		disp('Mult in function det3x3 10 overflow');
        disp({Mult_DetM_3x3_n_10_int(i)});
    end
	%% Проверка выходной разрядности умножителя
    if (Mult_DetM_3x3_n_7_total_width(i) > sim_options.width_hilbert-1)
		disp('Mult total width in function det3x3 10 higher than 64');
        disp({Mult_DetM_3x3_n_10_int(i)});
    end
end

%% 1 3x3
DetM_3x3_n_11 = Mult_DetM_2x2_n_10_a33 - Mult_DetM_2x2_n_9_a34 + Mult_DetM_2x2_n_8_a35; 
Det_3x3_LU_matlab(1) = (det(a(3:end,3:end))); % 3 4 5

DetM_3x3_n_12 = Mult_DetM_2x2_n_10_a32 - Mult_DetM_2x2_n_7_a34 + Mult_DetM_2x2_n_6_a35; 
Det_3x3_LU_matlab(2) = (det([a(3:end,2), a(3:end,4:5)])); % 2 4 5 

DetM_3x3_n_13 = Mult_DetM_2x2_n_9_a32 - Mult_DetM_2x2_n_7_a33 + Mult_DetM_2x2_n_5_a35;
Det_3x3_LU_matlab(3) = (det([a(3:end,2), a(3:end,3), a(3:end,5)])); % 2 3 5

DetM_3x3_n_14 = Mult_DetM_2x2_n_8_a32 - Mult_DetM_2x2_n_6_a33 + Mult_DetM_2x2_n_5_a34;
Det_3x3_LU_matlab(4) = (det([a(3:end,2), a(3:end,3), a(3:end,4)])); 

% 2 3x3
DetM_3x3_n_22 = Mult_DetM_2x2_n_10_a31 - Mult_DetM_2x2_n_4_a34 + Mult_DetM_2x2_n_3_a35; 
Det_3x3_LU_matlab(5) = (det([a(3:end,1), a(3:end,4), a(3:end,5)])); 

DetM_3x3_n_23 = Mult_DetM_2x2_n_9_a31 - Mult_DetM_2x2_n_4_a33 + Mult_DetM_2x2_n_2_a35;
Det_3x3_LU_matlab(6) = (det([a(3:end,1), a(3:end,3), a(3:end,5)])); 

DetM_3x3_n_24 = Mult_DetM_2x2_n_8_a31 - Mult_DetM_2x2_n_3_a33 + Mult_DetM_2x2_n_2_a34;
Det_3x3_LU_matlab(8) = (det([a(3:end,1), a(3:end,2), a(3:end,5)])); 

% 3 3x3
DetM_3x3_n_33 = Mult_DetM_2x2_n_7_a31 - Mult_DetM_2x2_n_4_a32 + Mult_DetM_2x2_n_1_a35;
Det_3x3_LU_matlab(9) = (det([a(3:end,1), a(3:end,2), a(3:end,4)])); 

DetM_3x3_n_34 = Mult_DetM_2x2_n_6_a31 - Mult_DetM_2x2_n_3_a32 + Mult_DetM_2x2_n_1_a34;
Det_3x3_LU_matlab(7) = (det([a(3:end,1), a(3:end,3), a(3:end,4)])); 

% 4 3x3
DetM_3x3_n_44 = Mult_DetM_2x2_n_5_a31 - Mult_DetM_2x2_n_2_a32 + Mult_DetM_2x2_n_1_a33;
Det_3x3_LU_matlab(10) = (det([a(3:end,1), a(3:end,2), a(3:end,3)])); 


%%
[DetM_3x3_n_44_int_sum1, DetM_3x3_n_44_int_sum1_overflow, DetM_3x3_n_44_int_sum1_abs, DetM_3x3_n_44_int_sum1_width_total] = ...
    adder(Mult_DetM_3x3_n_1_int(1), -Mult_DetM_3x3_n_2_int(1), sim_options.type_3x3_det, sim_options.width_hilbert);

	% Наложение маски на первый умножитель
    if sim_options.enable_mask == true
       c1 = bitmask(DetM_3x3_n_44_int_sum1, sim_options.type_3x3_det, width_mult1(i));
       if c1 ~= DetM_3x3_n_44_int_sum1
           disp('Bit mask error mult 1 det 2x2');
           % c2 = bitmask(sum(i+1,n), sim_options.int_size, width_sum(i+1));
           disp(width_total_mult1(i));
           disp({c1, DetM_3x3_n_44_int_sum1});
           disp({sim_options.freq, sim_options.SNR});
       end
       DetM_3x3_n_44_int_sum1(i) = c1;
    end
	
	%% Проверка переполнения сумматора
    if (DetM_3x3_n_44_int_sum1_overflow == 1)
		disp('Sum det3x3 44 overflow');
        disp({DetM_3x3_n_44_int_sum1});
    end
	%% Проверка выходной разрядности сумматора
    if (DetM_3x3_n_44_int_sum1_width_total > sim_options.width_hilbert)
		disp('Sum total width in function det3x3 44 higher than 64');
        disp({DetM_3x3_n_44_int_sum1});
    end
	
[DetM_3x3_n_44_int, DetM_3x3_n_44_int_overflow, DetM_3x3_n_44_int_abs, DetM_3x3_n_44_int_width_total] = ...
    adder(Mult_DetM_3x3_n_5_int(1), DetM_3x3_n_44_int_sum1, sim_options.type_3x3_det, sim_options.width_hilbert);
	
	% Наложение маски на первый умножитель
    if sim_options.enable_mask == true
       c1 = bitmask(DetM_3x3_n_44_int, sim_options.type_3x3_det, width_mult1);
       if c1 ~= DetM_3x3_n_44_int_sum1
           disp('Bit mask error mult 1 det 2x2');
           % c2 = bitmask(sum(i+1,n), sim_options.int_size, width_sum(i+1));
           disp(width_total_mult1);
           disp({c1, DetM_3x3_n_44_int_sum1});
           disp({sim_options.freq, sim_options.SNR});
       end
       DetM_3x3_n_44_int_sum1 = c1;
    end
	
	
    %% Проверка переполнения сумматора
    if (DetM_3x3_n_44_int_overflow == 1)
		disp('Sum det3x3 44 overflow');
        disp({DetM_3x3_n_44_int});
    end
	%% Проверка выходной разрядности сумматора
    if (DetM_3x3_n_44_int_width_total > sim_options.width_hilbert)
		disp('Sum total width in function det3x3 44 higher than 64');
        disp({DetM_3x3_n_44_int});
    end
    %%
[DetM_3x3_n_34_int_sum1, DetM_3x3_n_34_int_sum1_overflow, DetM_3x3_n_34_int_sum1_abs, DetM_3x3_n_34_int_sum1_width_total] = ...
    adder(Mult_DetM_3x3_n_1_int(2), -Mult_DetM_3x3_n_3_int(1), sim_options.type_3x3_det, sim_options.width_hilbert);
[DetM_3x3_n_34_int, DetM_3x3_n_34_int_overflow, DetM_3x3_n_34_int_abs, DetM_3x3_n_34_int_width_total] = ...
    adder(Mult_DetM_3x3_n_6_int(1), DetM_3x3_n_34_int_sum1, sim_options.type_3x3_det, sim_options.width_hilbert);
    %% Проверка переполнения сумматора
    if (DetM_3x3_n_34_int_overflow == 1)
		disp('Sum det3x3 34 overflow');
        disp({DetM_3x3_n_34_int});
    end
	%% Проверка выходной разрядности сумматора
    if (DetM_3x3_n_34_int_width_total > sim_options.width_hilbert)
		disp('Sum total width in function det3x3 34 higher than 64');
        disp({DetM_3x3_n_34_int});
    end
    %%
[DetM_3x3_n_33_int_sum1, DetM_3x3_n_33_int_sum1_overflow, DetM_3x3_n_33_int_sum1_abs, DetM_3x3_n_33_int_sum1_width_total] = ...
    adder(Mult_DetM_3x3_n_1_int(3), -Mult_DetM_3x3_n_4_int(1), sim_options.type_3x3_det, sim_options.width_hilbert);
[DetM_3x3_n_33_int, DetM_3x3_n_33_int_overflow, DetM_3x3_n_33_int_abs, DetM_3x3_n_33_int_width_total] = ...
    adder(Mult_DetM_3x3_n_7_int(1), DetM_3x3_n_33_int_sum1, sim_options.type_3x3_det, sim_options.width_hilbert);
    %% Проверка переполнения сумматора
    if (DetM_3x3_n_33_int_overflow == 1)
		disp('Sum det3x3 33 overflow');
        disp({DetM_3x3_n_33_int});
    end
	%% Проверка выходной разрядности сумматора
    if (DetM_3x3_n_33_int_width_total > sim_options.width_hilbert)
		disp('Sum total width in function det3x3 33 higher than 64');
        disp({DetM_3x3_n_33_int});
    end
    %%
[DetM_3x3_n_24_int_sum1, DetM_3x3_n_24_int_sum1_overflow, DetM_3x3_n_24_int_sum1_abs, DetM_3x3_n_24_int_sum1_width_total] = ...
    adder(Mult_DetM_3x3_n_2_int(2), -Mult_DetM_3x3_n_3_int(2), sim_options.type_3x3_det, sim_options.width_hilbert);
[DetM_3x3_n_24_int, DetM_3x3_n_24_int_overflow, DetM_3x3_n_24_int_abs, DetM_3x3_n_24_int_width_total] = ...
    adder(Mult_DetM_3x3_n_8_int(1), DetM_3x3_n_24_int_sum1, sim_options.type_3x3_det, sim_options.width_hilbert);
    %% Проверка переполнения сумматора
    if (DetM_3x3_n_24_int_overflow == 1)
		disp('Sum det3x3 24 overflow');
        disp({DetM_3x3_n_24_int});
    end
	%% Проверка выходной разрядности сумматора
    if (DetM_3x3_n_24_int_width_total > sim_options.width_hilbert)
		disp('Sum total width in function det3x3 24 higher than 64');
        disp({DetM_3x3_n_24_int});
    end
    %%
[DetM_3x3_n_23_int_sum1, DetM_3x3_n_23_int_sum1_overflow, DetM_3x3_n_23_int_sum1_abs, DetM_3x3_n_23_int_sum1_width_total] = ...
    adder(Mult_DetM_3x3_n_2_int(3), -Mult_DetM_3x3_n_4_int(2), sim_options.type_3x3_det, sim_options.width_hilbert);
[DetM_3x3_n_23_int, DetM_3x3_n_23_int_overflow, DetM_3x3_n_23_int_abs, DetM_3x3_n_23_int_width_total] = ...
    adder(Mult_DetM_3x3_n_9_int(1), DetM_3x3_n_23_int_sum1, sim_options.type_3x3_det, sim_options.width_hilbert);
    %% Проверка переполнения сумматора
    if (DetM_3x3_n_23_int_overflow == 1)
		disp('Sum det3x3 23 overflow');
        disp({DetM_3x3_n_23_int});
    end
	%% Проверка выходной разрядности сумматора
    if (DetM_3x3_n_23_int_width_total > sim_options.width_hilbert)
		disp('Sum total width in function det3x3 23 higher than 64');
        disp({DetM_3x3_n_23_int});
    end
    %%
[DetM_3x3_n_22_int_sum1, DetM_3x3_n_22_int_sum1_overflow, DetM_3x3_n_22_int_sum1_abs, DetM_3x3_n_22_int_sum1_width_total] = ...
    adder(Mult_DetM_3x3_n_3_int(3), -Mult_DetM_3x3_n_4_int(3), sim_options.type_3x3_det, sim_options.width_hilbert);
[DetM_3x3_n_22_int, DetM_3x3_n_22_int_overflow, DetM_3x3_n_22_int_abs, DetM_3x3_n_22_int_width_total] = ...
    adder(Mult_DetM_3x3_n_10_int(1), DetM_3x3_n_22_int_sum1, sim_options.type_3x3_det, sim_options.width_hilbert);
    %% Проверка переполнения сумматора
    if (DetM_3x3_n_22_int_overflow == 1)
		disp('Sum det3x3 22 overflow');
        disp({DetM_3x3_n_22_int});
    end
	%% Проверка выходной разрядности сумматора
    if (DetM_3x3_n_22_int_width_total > sim_options.width_hilbert-1)
		disp('Sum total width in function det3x3 22 higher than 64');
        disp({DetM_3x3_n_22_int});
    end
    %%
[DetM_3x3_n_14_int_sum1, DetM_3x3_n_14_int_sum1_overflow, DetM_3x3_n_14_int_sum1_abs, DetM_3x3_n_14_int_sum1_width_total] = ...
    adder(Mult_DetM_3x3_n_5_int(2), -Mult_DetM_3x3_n_6_int(2), sim_options.type_3x3_det, sim_options.width_hilbert);
[DetM_3x3_n_14_int, DetM_3x3_n_14_int_overflow, DetM_3x3_n_14_int_abs, DetM_3x3_n_14_int_width_total] = ...
    adder(Mult_DetM_3x3_n_8_int(2), DetM_3x3_n_14_int_sum1, sim_options.type_3x3_det, sim_options.width_hilbert);
    %% Проверка переполнения сумматора
    if (DetM_3x3_n_14_int_overflow == 1)
		disp('Sum det3x3 14 overflow');
        disp({DetM_3x3_n_14_int});
    end
	%% Проверка выходной разрядности сумматора
    if (DetM_3x3_n_14_int_width_total > sim_options.width_hilbert-1)
		disp('Sum total width in function det3x3 14 higher than 64');
        disp({DetM_3x3_n_14_int});
    end
    %%
[DetM_3x3_n_13_int_sum1, DetM_3x3_n_13_int_sum1_overflow, DetM_3x3_n_13_int_sum1_abs, DetM_3x3_n_13_int_sum1_width_total] = ...
    adder(Mult_DetM_3x3_n_5_int(3), -Mult_DetM_3x3_n_7_int(2), sim_options.type_3x3_det, sim_options.width_hilbert);
[DetM_3x3_n_13_int, DetM_3x3_n_13_int_overflow, DetM_3x3_n_13_int_abs, DetM_3x3_n_13_int_width_total] = ...
    adder(Mult_DetM_3x3_n_9_int(2), DetM_3x3_n_13_int_sum1, sim_options.type_3x3_det, sim_options.width_hilbert);
    %% Проверка переполнения сумматора
    if (DetM_3x3_n_13_int_overflow == 1)
		disp('Sum det3x3 13 overflow');
        disp({DetM_3x3_n_13_int});
    end
	%% Проверка выходной разрядности сумматора
    if (DetM_3x3_n_13_int_width_total > sim_options.width_hilbert-1)
		disp('Sum total width in function det3x3 13 higher than 64');
        disp({DetM_3x3_n_13_int});
    end
    %%
[DetM_3x3_n_12_int_sum1, DetM_3x3_n_12_int_sum1_overflow, DetM_3x3_n_12_int_sum1_abs, DetM_3x3_n_12_int_sum1_width_total] = ...
    adder(Mult_DetM_3x3_n_6_int(3), -Mult_DetM_3x3_n_7_int(3), sim_options.type_3x3_det, sim_options.width_hilbert);
[DetM_3x3_n_12_int, DetM_3x3_n_12_int_overflow, DetM_3x3_n_12_int_abs, DetM_3x3_n_12_int_width_total] = ...
    adder(Mult_DetM_3x3_n_10_int(2), DetM_3x3_n_12_int_sum1, sim_options.type_3x3_det, sim_options.width_hilbert);
    %% Проверка переполнения сумматора
    if (DetM_3x3_n_12_int_overflow == 1)
		disp('Sum det3x3 12 overflow');
        disp({DetM_3x3_n_12_int});
    end
	%% Проверка выходной разрядности сумматора
    if (DetM_3x3_n_12_int_width_total > sim_options.width_hilbert-1)
		disp('Sum total width in function det3x3 12 higher than 64');
        disp({DetM_3x3_n_12_int});
    end
    %%
[DetM_3x3_n_11_int_sum1, DetM_3x3_n_11_int_sum1_overflow, DetM_3x3_n_11_int_sum1_abs, DetM_3x3_n_11_int_sum1_width_total] = ...
    adder(Mult_DetM_3x3_n_8_int(3), -Mult_DetM_3x3_n_9_int(3), sim_options.type_3x3_det, sim_options.width_hilbert);
[DetM_3x3_n_11_int, DetM_3x3_n_11_int_overflow, DetM_3x3_n_11_int_abs, DetM_3x3_n_11_int_width_total] = ...
    adder(Mult_DetM_3x3_n_10_int(3), DetM_3x3_n_11_int_sum1, sim_options.type_3x3_det, sim_options.width_hilbert);
    %% Проверка переполнения сумматора
    if (DetM_3x3_n_11_int_overflow == 1)
		disp('Sum det3x3 11 overflow');
        disp({DetM_3x3_n_11_int});
    end
	%% Проверка выходной разрядности сумматора
    if (DetM_3x3_n_11_int_width_total > sim_options.width_hilbert-1)
		disp('Sum total width in function det3x3 11 higher than 64');
        disp({DetM_3x3_n_11_int});
    end
    %%
    if DetM_3x3_n_11 < 0
        DetM_3x3_array(1) = DetM_3x3_n_11 * -1;
    else
        DetM_3x3_array(1) = DetM_3x3_n_11;
    end

    if DetM_3x3_n_12 < 0
        DetM_3x3_array(2) = DetM_3x3_n_12 * -1;
    else
        DetM_3x3_array(2) = DetM_3x3_n_12;
    end

    if DetM_3x3_n_13 < 0
        DetM_3x3_array(3) = DetM_3x3_n_13 * -1;
    else
        DetM_3x3_array(3) = DetM_3x3_n_13;
    end
    if DetM_3x3_n_14 < 0
        DetM_3x3_array(4) = DetM_3x3_n_14 * -1;
    else
        DetM_3x3_array(4) = DetM_3x3_n_14;
    end
    if DetM_3x3_n_22 < 0
        DetM_3x3_array(5) = DetM_3x3_n_22 * -1;
    else
        DetM_3x3_array(5) = DetM_3x3_n_22;
    end
    if DetM_3x3_n_23 < 0
        DetM_3x3_array(6) = DetM_3x3_n_23 * -1;
    else
        DetM_3x3_array(6) = DetM_3x3_n_23;
    end
    if DetM_3x3_n_24 < 0
        DetM_3x3_array(7) = DetM_3x3_n_24 * -1;
    else
        DetM_3x3_array(7) = DetM_3x3_n_24;
    end
    if DetM_3x3_n_33 < 0
        DetM_3x3_array(8) = DetM_3x3_n_33 * -1;
    else
        DetM_3x3_array(8) = DetM_3x3_n_33;
    end
    if DetM_3x3_n_34 < 0
        DetM_3x3_array(9) = DetM_3x3_n_34 * -1;
    else
        DetM_3x3_array(9) = DetM_3x3_n_34;
    end
    if DetM_3x3_n_44 < 0
        DetM_3x3_array(10) = DetM_3x3_n_44 * -1;
    else
        DetM_3x3_array(10) = DetM_3x3_n_44;
    end
%%

DetM_3x3_int_pre_sum_array(1) = DetM_3x3_n_11_int_sum1;
DetM_3x3_int_pre_sum_array(2) = DetM_3x3_n_12_int_sum1;
DetM_3x3_int_pre_sum_array(3) = DetM_3x3_n_13_int_sum1;
DetM_3x3_int_pre_sum_array(4) = DetM_3x3_n_14_int_sum1;
DetM_3x3_int_pre_sum_array(5) = DetM_3x3_n_22_int_sum1;
DetM_3x3_int_pre_sum_array(6) = DetM_3x3_n_23_int_sum1;
DetM_3x3_int_pre_sum_array(7) = DetM_3x3_n_24_int_sum1;
DetM_3x3_int_pre_sum_array(8) = DetM_3x3_n_33_int_sum1;
DetM_3x3_int_pre_sum_array(9) = DetM_3x3_n_34_int_sum1;
DetM_3x3_int_pre_sum_array(10) = DetM_3x3_n_44_int_sum1;

DetM_3x3_int_pre_sum_width_total(1) = DetM_3x3_n_11_int_sum1_width_total;
DetM_3x3_int_pre_sum_width_total(2) = DetM_3x3_n_12_int_sum1_width_total;
DetM_3x3_int_pre_sum_width_total(3) = DetM_3x3_n_13_int_sum1_width_total;
DetM_3x3_int_pre_sum_width_total(4) = DetM_3x3_n_14_int_sum1_width_total;
DetM_3x3_int_pre_sum_width_total(5) = DetM_3x3_n_22_int_sum1_width_total;
DetM_3x3_int_pre_sum_width_total(6) = DetM_3x3_n_23_int_sum1_width_total;
DetM_3x3_int_pre_sum_width_total(7) = DetM_3x3_n_24_int_sum1_width_total;
DetM_3x3_int_pre_sum_width_total(8) = DetM_3x3_n_33_int_sum1_width_total;
DetM_3x3_int_pre_sum_width_total(9) = DetM_3x3_n_34_int_sum1_width_total;
DetM_3x3_int_pre_sum_width_total(10) = DetM_3x3_n_44_int_sum1_width_total;

DetM_3x3_int_sum_array(1) = DetM_3x3_n_11_int_abs;
DetM_3x3_int_sum_array(2) = DetM_3x3_n_12_int_abs;
DetM_3x3_int_sum_array(3) = DetM_3x3_n_13_int_abs;
DetM_3x3_int_sum_array(4) = DetM_3x3_n_14_int_abs;
DetM_3x3_int_sum_array(5) = DetM_3x3_n_22_int_abs;
DetM_3x3_int_sum_array(6) = DetM_3x3_n_23_int_abs;
DetM_3x3_int_sum_array(7) = DetM_3x3_n_24_int_abs;
DetM_3x3_int_sum_array(8) = DetM_3x3_n_33_int_abs;
DetM_3x3_int_sum_array(9) = DetM_3x3_n_34_int_abs;
DetM_3x3_int_sum_array(10) = DetM_3x3_n_44_int_abs;

DetM_3x3_int_sum_width_total(1) = DetM_3x3_n_11_int_width_total;
DetM_3x3_int_sum_width_total(2) = DetM_3x3_n_12_int_width_total;
DetM_3x3_int_sum_width_total(3) = DetM_3x3_n_13_int_width_total;
DetM_3x3_int_sum_width_total(4) = DetM_3x3_n_14_int_width_total;
DetM_3x3_int_sum_width_total(5) = DetM_3x3_n_22_int_width_total;
DetM_3x3_int_sum_width_total(6) = DetM_3x3_n_23_int_width_total;
DetM_3x3_int_sum_width_total(7) = DetM_3x3_n_24_int_width_total;
DetM_3x3_int_sum_width_total(8) = DetM_3x3_n_33_int_width_total;
DetM_3x3_int_sum_width_total(9) = DetM_3x3_n_34_int_width_total;
DetM_3x3_int_sum_width_total(10) = DetM_3x3_n_44_int_width_total;

                                                                    %% 4x4
% double                                                            
var_mult_a_det11 = a(2,2) * DetM_3x3_n_11; 
var_mult_a_det12 = a(2,3) * DetM_3x3_n_12;
var_mult_a_det13 = a(2,4) * DetM_3x3_n_13;
var_mult_a_det14 = a(2,5) * DetM_3x3_n_14;

DetM_4x4_n_1 = var_mult_a_det11 - var_mult_a_det12  + var_mult_a_det13 - var_mult_a_det14;
Det_4x4_LU_matlab(1) = det([a(2:end,2), a(2:end,3), a(2:end,4), a(2:end,5)]); 

var_mult_a_det21 = a(2,1) * DetM_3x3_n_11;
var_mult_a_det22 = a(2,3) * DetM_3x3_n_22;
var_mult_a_det23 = a(2,4) * DetM_3x3_n_23;
var_mult_a_det24 = a(2,5) * DetM_3x3_n_24;

DetM_4x4_n_2 = var_mult_a_det21 - var_mult_a_det22 + var_mult_a_det23 - var_mult_a_det24;
Det_4x4_LU_matlab(2) = det([a(2:end,1), a(2:end,3), a(2:end,4), a(2:end,5)]); 

var_mult_a_det31 = a(2,1) * DetM_3x3_n_12;
var_mult_a_det32 = a(2,2) * DetM_3x3_n_22;
var_mult_a_det33 = a(2,4) * DetM_3x3_n_33;
var_mult_a_det34 = a(2,5) * DetM_3x3_n_34;

DetM_4x4_n_3 = var_mult_a_det31 - var_mult_a_det32 + var_mult_a_det33 - var_mult_a_det34;
Det_4x4_LU_matlab(3) = det([a(2:end,1), a(2:end,2), a(2:end,4), a(2:end,5)]); 

var_mult_a_det41 = a(2,1) * DetM_3x3_n_13;
var_mult_a_det42 = a(2,2) * DetM_3x3_n_23;
var_mult_a_det43 = a(2,3) * DetM_3x3_n_33;
var_mult_a_det44 = a(2,5) * DetM_3x3_n_44;

DetM_4x4_n_4 = var_mult_a_det41 - var_mult_a_det42 + var_mult_a_det43 - var_mult_a_det44;
Det_4x4_LU_matlab(4) = det([a(2:end,1), a(2:end,2), a(2:end,3), a(2:end,5)]); 

var_mult_a_det51 = a(2,1) * DetM_3x3_n_14;
var_mult_a_det52 = a(2,2) * DetM_3x3_n_24;
var_mult_a_det53 = a(2,3) * DetM_3x3_n_34;
var_mult_a_det54 = a(2,4) * DetM_3x3_n_44;

DetM_4x4_n_5 = var_mult_a_det51 - var_mult_a_det52 + var_mult_a_det53 - var_mult_a_det54;
Det_4x4_LU_matlab(5) = det([a(2:end,1), a(2:end,2), a(2:end,3), a(2:end,4)]); 

DetM_4x4_array(1) = (DetM_4x4_n_1);
DetM_4x4_array(2) = (DetM_4x4_n_2);
DetM_4x4_array(3) = (DetM_4x4_n_3);
DetM_4x4_array(4) = (DetM_4x4_n_4);
DetM_4x4_array(5) = (DetM_4x4_n_5);

for n = 1:5
    if Det_4x4_LU_matlab(n) < 0
        Det_4x4_LU_matlab(n) = Det_4x4_LU_matlab(n) * -1;
    end
    if DetM_4x4_array(n) < 0
        DetM_4x4_array(n) = DetM_4x4_array(n) * -1;
    end
end

%% integer
[var_mult_a_det11_int, var_mult_a_det11_int_overflow, var_mult_a_det11_abs, var_mult_a_det11_width_total] = mult(cast(a_int(2,2), sim_options.type_4x4_det), cast(DetM_3x3_n_11_int, sim_options.type_4x4_det), sim_options.type_4x4_det, sim_options.width_hilbert);
[var_mult_a_det12_int, var_mult_a_det12_int_overflow, var_mult_a_det12_abs, var_mult_a_det12_width_total] = mult(cast(a_int(2,3), sim_options.type_4x4_det), cast(DetM_3x3_n_12_int, sim_options.type_4x4_det), sim_options.type_4x4_det, sim_options.width_hilbert);
[var_mult_a_det13_int, var_mult_a_det13_int_overflow, var_mult_a_det13_abs, var_mult_a_det13_width_total] = mult(cast(a_int(2,4), sim_options.type_4x4_det), cast(DetM_3x3_n_13_int, sim_options.type_4x4_det), sim_options.type_4x4_det, sim_options.width_hilbert);
[var_mult_a_det14_int, var_mult_a_det14_int_overflow, var_mult_a_det14_abs, var_mult_a_det14_width_total] = mult(cast(a_int(2,5), sim_options.type_4x4_det), cast(DetM_3x3_n_14_int, sim_options.type_4x4_det), sim_options.type_4x4_det, sim_options.width_hilbert);
[var_mult_a_det14_int, var_mult_a_det14_int_overflow, var_mult_a_det14_abs, var_mult_a_det14_width_total] = mult(cast(a_int(2,5), sim_options.type_4x4_det), cast(DetM_3x3_n_14_int, sim_options.type_4x4_det), sim_options.type_4x4_det, sim_options.width_hilbert);

[var_mult_a_det21_int, var_mult_a_det21_int_overflow, var_mult_a_det21_abs, var_mult_a_det21_width_total] = mult(cast(a_int(2,1), sim_options.type_4x4_det), cast(DetM_3x3_n_11_int, sim_options.type_4x4_det), sim_options.type_4x4_det, sim_options.width_hilbert);
[var_mult_a_det22_int, var_mult_a_det22_int_overflow, var_mult_a_det22_abs, var_mult_a_det22_width_total] = mult(cast(a_int(2,3), sim_options.type_4x4_det), cast(DetM_3x3_n_22_int, sim_options.type_4x4_det), sim_options.type_4x4_det, sim_options.width_hilbert);
[var_mult_a_det23_int, var_mult_a_det23_int_overflow, var_mult_a_det23_abs, var_mult_a_det23_width_total] = mult(cast(a_int(2,4), sim_options.type_4x4_det), cast(DetM_3x3_n_23_int, sim_options.type_4x4_det), sim_options.type_4x4_det, sim_options.width_hilbert);
[var_mult_a_det24_int, var_mult_a_det24_int_overflow, var_mult_a_det24_abs, var_mult_a_det24_width_total] = mult(cast(a_int(2,5), sim_options.type_4x4_det), cast(DetM_3x3_n_24_int, sim_options.type_4x4_det), sim_options.type_4x4_det, sim_options.width_hilbert);

[var_mult_a_det31_int, var_mult_a_det31_int_overflow, var_mult_a_det31_abs, var_mult_a_det31_width_total] = mult(cast(a_int(2,1), sim_options.type_4x4_det), cast(DetM_3x3_n_12_int, sim_options.type_4x4_det), sim_options.type_4x4_det, sim_options.width_hilbert);
[var_mult_a_det32_int, var_mult_a_det32_int_overflow, var_mult_a_det32_abs, var_mult_a_det32_width_total] = mult(cast(a_int(2,2), sim_options.type_4x4_det), cast(DetM_3x3_n_22_int, sim_options.type_4x4_det), sim_options.type_4x4_det, sim_options.width_hilbert);
[var_mult_a_det33_int, var_mult_a_det33_int_overflow, var_mult_a_det33_abs, var_mult_a_det33_width_total] = mult(cast(a_int(2,4), sim_options.type_4x4_det), cast(DetM_3x3_n_33_int, sim_options.type_4x4_det), sim_options.type_4x4_det, sim_options.width_hilbert);
[var_mult_a_det34_int, var_mult_a_det34_int_overflow, var_mult_a_det34_abs, var_mult_a_det34_width_total] = mult(cast(a_int(2,5), sim_options.type_4x4_det), cast(DetM_3x3_n_34_int, sim_options.type_4x4_det), sim_options.type_4x4_det, sim_options.width_hilbert);

[var_mult_a_det41_int, var_mult_a_det41_int_overflow, var_mult_a_det41_abs, var_mult_a_det41_width_total] = mult(cast(a_int(2,1), sim_options.type_4x4_det), cast(DetM_3x3_n_13_int, sim_options.type_4x4_det), sim_options.type_4x4_det, sim_options.width_hilbert);
[var_mult_a_det42_int, var_mult_a_det42_int_overflow, var_mult_a_det42_abs, var_mult_a_det42_width_total] = mult(cast(a_int(2,2), sim_options.type_4x4_det), cast(DetM_3x3_n_23_int, sim_options.type_4x4_det), sim_options.type_4x4_det, sim_options.width_hilbert);
[var_mult_a_det43_int, var_mult_a_det43_int_overflow, var_mult_a_det43_abs, var_mult_a_det43_width_total] = mult(cast(a_int(2,3), sim_options.type_4x4_det), cast(DetM_3x3_n_33_int, sim_options.type_4x4_det), sim_options.type_4x4_det, sim_options.width_hilbert);
[var_mult_a_det44_int, var_mult_a_det44_int_overflow, var_mult_a_det44_abs, var_mult_a_det44_width_total] = mult(cast(a_int(2,5), sim_options.type_4x4_det), cast(DetM_3x3_n_44_int, sim_options.type_4x4_det), sim_options.type_4x4_det, sim_options.width_hilbert);

[var_mult_a_det51_int, var_mult_a_det51_int_overflow, var_mult_a_det51_abs, var_mult_a_det51_width_total] = mult(cast(a_int(2,1), sim_options.type_4x4_det), cast(DetM_3x3_n_14_int, sim_options.type_4x4_det), sim_options.type_4x4_det, sim_options.width_hilbert);
[var_mult_a_det52_int, var_mult_a_det52_int_overflow, var_mult_a_det52_abs, var_mult_a_det52_width_total] = mult(cast(a_int(2,2), sim_options.type_4x4_det), cast(DetM_3x3_n_24_int, sim_options.type_4x4_det), sim_options.type_4x4_det, sim_options.width_hilbert);
[var_mult_a_det53_int, var_mult_a_det53_int_overflow, var_mult_a_det53_abs, var_mult_a_det53_width_total] = mult(cast(a_int(2,3), sim_options.type_4x4_det), cast(DetM_3x3_n_34_int, sim_options.type_4x4_det), sim_options.type_4x4_det, sim_options.width_hilbert);
[var_mult_a_det54_int, var_mult_a_det54_int_overflow, var_mult_a_det54_abs, var_mult_a_det54_width_total] = mult(cast(a_int(2,4), sim_options.type_4x4_det), cast(DetM_3x3_n_44_int, sim_options.type_4x4_det), sim_options.type_4x4_det, sim_options.width_hilbert);

%%
[DetM_4x4_n_1_int_sum1, DetM_4x4_n_1_int_sum1_overflow, DetM_4x4_n_1_int_sum1_abs, DetM_4x4_n_1_int_sum1_width_total] = adder(var_mult_a_det11_int,  -var_mult_a_det12_int,  sim_options.type_4x4_det, sim_options.width_hilbert);
[DetM_4x4_n_1_int_sum2, DetM_4x4_n_1_int_sum2_overflow, DetM_4x4_n_1_int_sum2_abs, DetM_4x4_n_1_int_sum2_width_total] = adder(var_mult_a_det13_int,  -var_mult_a_det14_int,  sim_options.type_4x4_det, sim_options.width_hilbert);
[DetM_4x4_n_1_int, 		DetM_4x4_n_1_int_overflow, 		DetM_4x4_n_1_int_abs, 		DetM_4x4_n_1_int_width_total] =		adder(DetM_4x4_n_1_int_sum1,  DetM_4x4_n_1_int_sum2, sim_options.type_4x4_det, sim_options.width_hilbert);

[DetM_4x4_n_2_int_sum1, DetM_4x4_n_2_int_sum1_overflow, DetM_4x4_n_2_int_sum1_abs, DetM_4x4_n_2_int_sum1_width_total] = adder(var_mult_a_det21_int,  -var_mult_a_det22_int,  sim_options.type_4x4_det, sim_options.width_hilbert);
[DetM_4x4_n_2_int_sum2, DetM_4x4_n_2_int_sum2_overflow, DetM_4x4_n_2_int_sum2_abs, DetM_4x4_n_2_int_sum2_width_total] = adder(var_mult_a_det23_int,  -var_mult_a_det24_int,  sim_options.type_4x4_det, sim_options.width_hilbert);
[DetM_4x4_n_2_int, 		DetM_4x4_n_2_int_overflow,		DetM_4x4_n_2_int_abs, 		DetM_4x4_n_2_int_width_total] =		adder(DetM_4x4_n_2_int_sum1,  DetM_4x4_n_2_int_sum2, sim_options.type_4x4_det, sim_options.width_hilbert);

[DetM_4x4_n_3_int_sum1, DetM_4x4_n_3_int_sum1_overflow, DetM_4x4_n_3_int_sum1_abs, DetM_4x4_n_3_int_sum1_width_total] = adder(var_mult_a_det31_int,  -var_mult_a_det32_int,  sim_options.type_4x4_det, sim_options.width_hilbert);
[DetM_4x4_n_3_int_sum2, DetM_4x4_n_3_int_sum2_overflow, DetM_4x4_n_3_int_sum2_abs, DetM_4x4_n_3_int_sum2_width_total] = adder(var_mult_a_det33_int,  -var_mult_a_det34_int,  sim_options.type_4x4_det, sim_options.width_hilbert);
[DetM_4x4_n_3_int, 		DetM_4x4_n_3_int_overflow,		DetM_4x4_n_3_int_abs, 		DetM_4x4_n_3_int_width_total] = 	adder(DetM_4x4_n_3_int_sum1,  DetM_4x4_n_3_int_sum2, sim_options.type_4x4_det, sim_options.width_hilbert);

[DetM_4x4_n_4_int_sum1, DetM_4x4_n_4_int_sum1_overflow, DetM_4x4_n_4_int_sum1_abs, DetM_4x4_n_4_int_sum1_width_total] = adder(var_mult_a_det41_int,  -var_mult_a_det42_int,  sim_options.type_4x4_det, sim_options.width_hilbert);
[DetM_4x4_n_4_int_sum2, DetM_4x4_n_4_int_sum2_overflow, DetM_4x4_n_4_int_sum2_abs, DetM_4x4_n_4_int_sum2_width_total] = adder(var_mult_a_det43_int,  -var_mult_a_det44_int,  sim_options.type_4x4_det, sim_options.width_hilbert);
[DetM_4x4_n_4_int, 		DetM_4x4_n_4_int_overflow,		DetM_4x4_n_4_int_abs, 		DetM_4x4_n_4_int_width_total] = 	adder(DetM_4x4_n_4_int_sum1,  DetM_4x4_n_4_int_sum2, sim_options.type_4x4_det, sim_options.width_hilbert);

[DetM_4x4_n_5_int_sum1, DetM_4x4_n_5_int_sum1_overflow, DetM_4x4_n_5_int_sum1_abs, DetM_4x4_n_5_int_sum1_width_total] = adder(var_mult_a_det51_int,  -var_mult_a_det52_int,  sim_options.type_4x4_det, sim_options.width_hilbert);
[DetM_4x4_n_5_int_sum2, DetM_4x4_n_5_int_sum2_overflow, DetM_4x4_n_5_int_sum2_abs, DetM_4x4_n_5_int_sum2_width_total] = adder(var_mult_a_det53_int,  -var_mult_a_det54_int,  sim_options.type_4x4_det, sim_options.width_hilbert);
[DetM_4x4_n_5_int, 		DetM_4x4_n_5_int_overflow,		DetM_4x4_n_5_int_abs, 		DetM_4x4_n_5_int_width_total] = 	adder(DetM_4x4_n_5_int_sum1,  DetM_4x4_n_5_int_sum2, sim_options.type_4x4_det, sim_options.width_hilbert);


DetM_4x4_int_mult_array(1) = var_mult_a_det11_abs;
DetM_4x4_int_mult_array(2) = var_mult_a_det12_abs;
DetM_4x4_int_mult_array(3) = var_mult_a_det13_abs;
DetM_4x4_int_mult_array(4) = var_mult_a_det14_abs;
DetM_4x4_int_mult_array(5) = var_mult_a_det21_abs;
DetM_4x4_int_mult_array(6) = var_mult_a_det22_abs;
DetM_4x4_int_mult_array(7) = var_mult_a_det23_abs;
DetM_4x4_int_mult_array(8) = var_mult_a_det24_abs;
DetM_4x4_int_mult_array(9) = var_mult_a_det31_abs;
DetM_4x4_int_mult_array(10) = var_mult_a_det32_abs;
DetM_4x4_int_mult_array(11) = var_mult_a_det33_abs;
DetM_4x4_int_mult_array(12) = var_mult_a_det34_abs;
DetM_4x4_int_mult_array(13) = var_mult_a_det41_abs;
DetM_4x4_int_mult_array(14) = var_mult_a_det42_abs;
DetM_4x4_int_mult_array(15) = var_mult_a_det43_abs;
DetM_4x4_int_mult_array(16) = var_mult_a_det44_abs;
DetM_4x4_int_mult_array(17) = var_mult_a_det51_abs;
DetM_4x4_int_mult_array(18) = var_mult_a_det52_abs;
DetM_4x4_int_mult_array(19) = var_mult_a_det53_abs;
DetM_4x4_int_mult_array(20) = var_mult_a_det54_abs;

DetM_4x4_int_mult_width_total(1) = var_mult_a_det11_width_total;
DetM_4x4_int_mult_width_total(2) = var_mult_a_det12_width_total;
DetM_4x4_int_mult_width_total(3) = var_mult_a_det13_width_total;
DetM_4x4_int_mult_width_total(4) = var_mult_a_det14_width_total;
DetM_4x4_int_mult_width_total(5) = var_mult_a_det21_width_total;
DetM_4x4_int_mult_width_total(6) = var_mult_a_det22_width_total;
DetM_4x4_int_mult_width_total(7) = var_mult_a_det23_width_total;
DetM_4x4_int_mult_width_total(8) = var_mult_a_det24_width_total;
DetM_4x4_int_mult_width_total(9) = var_mult_a_det31_width_total;
DetM_4x4_int_mult_width_total(10) = var_mult_a_det32_width_total;
DetM_4x4_int_mult_width_total(11) = var_mult_a_det33_width_total;
DetM_4x4_int_mult_width_total(12) = var_mult_a_det34_width_total;
DetM_4x4_int_mult_width_total(13) = var_mult_a_det41_width_total;
DetM_4x4_int_mult_width_total(14) = var_mult_a_det42_width_total;
DetM_4x4_int_mult_width_total(15) = var_mult_a_det43_width_total;
DetM_4x4_int_mult_width_total(16) = var_mult_a_det44_width_total;
DetM_4x4_int_mult_width_total(17) = var_mult_a_det51_width_total;
DetM_4x4_int_mult_width_total(18) = var_mult_a_det52_width_total;
DetM_4x4_int_mult_width_total(19) = var_mult_a_det53_width_total;
DetM_4x4_int_mult_width_total(20) = var_mult_a_det54_width_total;

DetM_4x4_int_pre_sum_array(1) = DetM_4x4_n_1_int_sum1;
DetM_4x4_int_pre_sum_array(2) = DetM_4x4_n_1_int_sum2;
DetM_4x4_int_pre_sum_array(3) = DetM_4x4_n_2_int_sum1;
DetM_4x4_int_pre_sum_array(4) = DetM_4x4_n_2_int_sum2;
DetM_4x4_int_pre_sum_array(5) = DetM_4x4_n_3_int_sum1;
DetM_4x4_int_pre_sum_array(6) = DetM_4x4_n_3_int_sum2;
DetM_4x4_int_pre_sum_array(7) = DetM_4x4_n_4_int_sum1;
DetM_4x4_int_pre_sum_array(8) = DetM_4x4_n_4_int_sum2;
DetM_4x4_int_pre_sum_array(9) = DetM_4x4_n_5_int_sum1;
DetM_4x4_int_pre_sum_array(10) = DetM_4x4_n_5_int_sum2;

DetM_4x4_int_pre_sum_width_total(1) = DetM_4x4_n_1_int_sum1_width_total;
DetM_4x4_int_pre_sum_width_total(2) = DetM_4x4_n_1_int_sum2_width_total;
DetM_4x4_int_pre_sum_width_total(3) = DetM_4x4_n_2_int_sum1_width_total;
DetM_4x4_int_pre_sum_width_total(4) = DetM_4x4_n_2_int_sum2_width_total;
DetM_4x4_int_pre_sum_width_total(5) = DetM_4x4_n_3_int_sum1_width_total;
DetM_4x4_int_pre_sum_width_total(6) = DetM_4x4_n_3_int_sum2_width_total;
DetM_4x4_int_pre_sum_width_total(7) = DetM_4x4_n_4_int_sum1_width_total;
DetM_4x4_int_pre_sum_width_total(8) = DetM_4x4_n_4_int_sum2_width_total;
DetM_4x4_int_pre_sum_width_total(9) = DetM_4x4_n_5_int_sum1_width_total;
DetM_4x4_int_pre_sum_width_total(10) = DetM_4x4_n_5_int_sum2_width_total;

DetM_4x4_int_sum_array_width_total(1) = DetM_4x4_n_1_int_width_total;
DetM_4x4_int_sum_array_width_total(2) = DetM_4x4_n_2_int_width_total;
DetM_4x4_int_sum_array_width_total(3) = DetM_4x4_n_3_int_width_total;
DetM_4x4_int_sum_array_width_total(4) = DetM_4x4_n_4_int_width_total;
DetM_4x4_int_sum_array_width_total(5) = DetM_4x4_n_5_int_width_total;

DetM_4x4_int_sum_array(1) = DetM_4x4_n_1_int_abs;
DetM_4x4_int_sum_array(2) = DetM_4x4_n_2_int_abs;
DetM_4x4_int_sum_array(3) = DetM_4x4_n_3_int_abs;
DetM_4x4_int_sum_array(4) = DetM_4x4_n_4_int_abs;
DetM_4x4_int_sum_array(5) = DetM_4x4_n_5_int_abs;

                                                                    %% 5x5
% double                                                                    
var1 = a(1,1) * DetM_4x4_n_1;
var2 = a(1,2) * DetM_4x4_n_2;
var3 = a(1,3) * DetM_4x4_n_3;
var4 = a(1,4) * DetM_4x4_n_4;
var5 = a(1,5) * DetM_4x4_n_5;
DetM_5x5 = var1 - var2 + var3 - var4 + var5;

%% integer
[var1_int, var1_int_overflow, var1_int_abs, var1_int_mult_width_total] = mult(double(a_int(1,1)), DetM_4x4_n_1_int, sim_options.type_5x5_det, sim_options.width_hilbert);
[var2_int, var2_int_overflow, var2_int_abs, var2_int_mult_width_total] = mult(double(a_int(1,2)), DetM_4x4_n_2_int, sim_options.type_5x5_det, sim_options.width_hilbert);
[var3_int, var3_int_overflow, var3_int_abs, var3_int_mult_width_total] = mult(double(a_int(1,3)), DetM_4x4_n_3_int, sim_options.type_5x5_det, sim_options.width_hilbert);
[var4_int, var4_int_overflow, var4_int_abs, var4_int_mult_width_total] = mult(double(a_int(1,4)), DetM_4x4_n_4_int, sim_options.type_5x5_det, sim_options.width_hilbert);
[var5_int, var5_int_overflow, var5_int_abs, var5_int_mult_width_total] = mult(double(a_int(1,5)), DetM_4x4_n_5_int, sim_options.type_5x5_det, sim_options.width_hilbert);

[DetM_5x5_int_sum1, DetM_5x5_int_sum1_overflow, DetM_5x5_int_sum1_abs, DetM_5x5_int_sum1_width_total] = adder(var1_int,  -var2_int, sim_options.type_5x5_det, sim_options.width_hilbert);
[DetM_5x5_int_sum2, DetM_5x5_int_sum2_overflow, DetM_5x5_int_sum2_abs, DetM_5x5_int_sum2_width_total] = adder(var3_int,  -var4_int, sim_options.type_5x5_det, sim_options.width_hilbert);
[DetM_5x5_int_sum3, DetM_5x5_int_sum3_overflow, DetM_5x5_int_sum3_abs, DetM_5x5_int_sum3_width_total] = adder(DetM_5x5_int_sum1,  DetM_5x5_int_sum2, 	sim_options.type_5x5_det, sim_options.width_hilbert);
[DetM_5x5_int, DetM_5x5_int_overflow, DetM_5x5_int_abs, DetM_5x5_int_width_total] = 					adder(DetM_5x5_int_sum3,  var5_int, 			sim_options.type_5x5_det, sim_options.width_hilbert);

DetM_5x5_int_mult_array(1) = var1_int_abs;
DetM_5x5_int_mult_array(2) = var2_int_abs;
DetM_5x5_int_mult_array(3) = var3_int_abs;
DetM_5x5_int_mult_array(4) = var4_int_abs;
DetM_5x5_int_mult_array(5) = var5_int_abs;

DetM_5x5_int_mult_array_width_total(1) = var1_int_mult_width_total;
DetM_5x5_int_mult_array_width_total(2) = var2_int_mult_width_total;
DetM_5x5_int_mult_array_width_total(3) = var3_int_mult_width_total;
DetM_5x5_int_mult_array_width_total(4) = var4_int_mult_width_total;
DetM_5x5_int_mult_array_width_total(5) = var5_int_mult_width_total;

DetM_5x5_int_pre_sum1_array(1) = DetM_5x5_int_sum1_abs;
DetM_5x5_int_pre_sum1_array(2) = DetM_5x5_int_sum2_abs;

DetM_5x5_int_pre_sum1_array_width_total(1) = DetM_5x5_int_sum1_width_total;
DetM_5x5_int_pre_sum1_array_width_total(2) = DetM_5x5_int_sum2_width_total;

%%

DetM_2x2_abs = DetM_2x2;
for n = 1:sim_options.num_det2x2
    if DetM_2x2_abs(n) < 0
        DetM_2x2_abs(n) = DetM_2x2_abs(n) * -1;
    end
end

% запись максимальных значений умножителей определителя 2х2
DetM_2x2_multiplier_total_abs(1:2:end) = Det2x2_mult1_abs;
DetM_2x2_multiplier_total_abs(2:2:end) = Det2x2_mult2_abs;

%%
Mult_DetM_3x3_array(1:3) = Mult_DetM_3x3_n_1_abs;
Mult_DetM_3x3_array(4:6) = Mult_DetM_3x3_n_2_abs;
Mult_DetM_3x3_array(7:9) = Mult_DetM_3x3_n_3_abs;
Mult_DetM_3x3_array(10:12) = Mult_DetM_3x3_n_4_abs;
Mult_DetM_3x3_array(13:15) = Mult_DetM_3x3_n_5_abs;
Mult_DetM_3x3_array(16:18) = Mult_DetM_3x3_n_6_abs;
Mult_DetM_3x3_array(19:21) = Mult_DetM_3x3_n_7_abs;
Mult_DetM_3x3_array(22:24) = Mult_DetM_3x3_n_8_abs;
Mult_DetM_3x3_array(25:27) = Mult_DetM_3x3_n_9_abs;
Mult_DetM_3x3_array(28:30) = Mult_DetM_3x3_n_10_abs;

Mult_DetM_3x3_array_mult_total_width(1:3) =   Mult_DetM_3x3_n_1_total_width;
Mult_DetM_3x3_array_mult_total_width(4:6) =   Mult_DetM_3x3_n_2_total_width;
Mult_DetM_3x3_array_mult_total_width(7:9) =   Mult_DetM_3x3_n_3_total_width;
Mult_DetM_3x3_array_mult_total_width(10:12) = Mult_DetM_3x3_n_4_total_width;
Mult_DetM_3x3_array_mult_total_width(13:15) = Mult_DetM_3x3_n_5_total_width;
Mult_DetM_3x3_array_mult_total_width(16:18) = Mult_DetM_3x3_n_6_total_width;
Mult_DetM_3x3_array_mult_total_width(19:21) = Mult_DetM_3x3_n_7_total_width;
Mult_DetM_3x3_array_mult_total_width(22:24) = Mult_DetM_3x3_n_8_total_width;
Mult_DetM_3x3_array_mult_total_width(25:27) = Mult_DetM_3x3_n_9_total_width;
Mult_DetM_3x3_array_mult_total_width(28:30) = Mult_DetM_3x3_n_10_total_width;

for n = 1:sim_options.num_det2x2
    if Det_3x3_LU_matlab(n) < 0
        Det_3x3_LU_matlab(n) = Det_3x3_LU_matlab(n) * -1;
    end
end
	

end