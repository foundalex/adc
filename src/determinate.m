
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
	DetM_5x5_int_width_total, ...
    s ...
] = determinate(data_in, data_in_int, sim_options)

DetM_2x2 = zeros(sim_options.num_det2x2,1);
DetM_2x2_int = cast(zeros(sim_options.num_det2x2,1), sim_options.type_2x2_det);

Det_2x2_LU_matlab = zeros(sim_options.num_det2x2,1);

Det2x2_mult1_abs = cast(zeros(sim_options.num_det2x2,1), sim_options.type_2x2_det);
Det2x2_mult2_abs = cast(zeros(sim_options.num_det2x2,1), sim_options.type_2x2_det);
s.Det2x2_sum_abs = cast(zeros(sim_options.num_det2x2,1), 	 sim_options.type_2x2_det);
s.DetM_2x2_multiplier_total_abs = cast(zeros(sim_options.num_det2x2*2,1), sim_options.type_2x2_det);

%% Определитель 3х3
Mult_DetM_3x3_int = cast(zeros(sim_options.num_det2x2*3,1), sim_options.type_3x3_det);

%%
Det_3x3_LU_matlab 						= 	zeros(sim_options.num_det2x2,1);
s.Mult_DetM_3x3_array 					= 	cast(zeros(sim_options.num_det2x2*3,1), sim_options.type_3x3_det);
s.Mult_DetM_3x3_array_mult_total_width 	= 	cast(zeros(sim_options.num_det2x2*3,1), sim_options.type_3x3_det);
s.DetM_3x3_int_pre_sum_array 				= 	cast(zeros(sim_options.num_det2x2,1), sim_options.type_3x3_det);
s.DetM_3x3_int_pre_sum_width_total 		= 	cast(zeros(sim_options.num_det2x2,1), sim_options.type_3x3_det);
s.DetM_3x3_int_sum_array 					= 	cast(zeros(sim_options.num_det2x2,1), sim_options.type_3x3_det);
s.DetM_3x3_int_sum_width_total 			=	cast(zeros(sim_options.num_det2x2,1), sim_options.type_3x3_det);

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

%%
a = data_in;
a_int = data_in_int;
                                                                            %% 2x2
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

%% Находим определители матриц 2x2
for i = 1:sim_options.num_det2x2
    % double LU

    ee = det(s.a{i});
    if (ee < 0)
        Det_2x2_LU_matlab(i) = ee * -1;
    else
        Det_2x2_LU_matlab(i) = ee;
    end

    %% Определитель 2х2

    DetM_2x2(i) = s.a{i}(1,1) * s.a{i}(2,2) - s.a{i}(2,1) * s.a{i}(1,2);

    DetM_2x2_abs(i) = DetM_2x2(i);
    if DetM_2x2_abs(i) < 0
        DetM_2x2_abs(i) = DetM_2x2_abs(i) * -1;
    end


    %% Умножители определителя 2х2
    [mult1_int(i), mult1_overflow(i), Det2x2_mult1_abs(i), width_total_mult1(i)] = mult(cast(s.a_int{i}(1,1),sim_options.type_2x2_det), ...
        cast(s.a_int{i}(2,2),sim_options.type_2x2_det), sim_options.type_2x2_det, sim_options.width_fractional);

    %% Проверка переполнения умножителя
    if (mult1_overflow(i) == 1)
	    disp('Mult1 in function det2x2 overflow');
        disp(sim_options.width_hilbert);
        disp({i});
    end
	% Проверка выходной разрядности умножителя
    if (width_total_mult1(i) > sim_options.width_hilbert)
	    disp('Mult1 total width in function det2x2 higher than 64');
        disp(sim_options.width_hilbert);
        disp({i});
    end

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
    end

    [mult2_int(i), mult2_overflow(i), Det2x2_mult2_abs(i), width_total_mult2(i)] = mult(cast(s.a_int{i}(2,1),sim_options.type_2x2_det), ...
        cast(s.a_int{i}(1,2),sim_options.type_2x2_det), sim_options.type_2x2_det, sim_options.width_fractional);

    %% Проверка переполнения умножителя
    if (mult2_overflow(i) == 1)
		disp('Mult2 in function det2x2 overflow');
        disp(sim_options.width_hilbert);
        disp({i});
    end
	% Проверка выходной разрядности умножителя
    if (width_total_mult2(i) > sim_options.width_hilbert)
		disp('Mult2 total width in function det2x2 higher than 64');
        disp(sim_options.width_hilbert);
        disp({i});
    end

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
    end

    %% Сумматоры определителя 2х2
	[DetM_2x2_int(i), sum_overflow(i), s.Det2x2_sum_abs{i}, width_total_sum(i)] = adder(cast(mult1_int(i),sim_options.type_2x2_det), ...
        cast(-mult2_int(i),sim_options.type_2x2_det), sim_options.type_2x2_det, sim_options.width_fractional);

    % Проверка переполнения сумматора
    if (sum_overflow(i) == 1)
		disp('Adder in function det2x2 overflow');
        disp(sim_options.width_hilbert);
        disp({s.a_int{i}});
    end
	% Проверка выходной разрядности сумматора
    if (width_total_sum(i) > sim_options.width_hilbert)
		disp('Sum total width in function det2x2 higher than 64');
        disp(sim_options.width_hilbert);
        disp({width_total_sum(i)});
    end
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
    end
end

% запись максимальных значений умножителей определителя 2х2
s.DetM_2x2_multiplier_total_abs(1:2:end) = Det2x2_mult1_abs;
s.DetM_2x2_multiplier_total_abs(2:2:end) = Det2x2_mult2_abs;

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

%% Умножители определителя 3х3

index(:,10) = [1, 2, 3];
index(:,9) = [1, 2, 4];
index(:,8) = [1, 2, 5];
index(:,7) = [1, 3, 4];
index(:,6) = [1, 3, 5];
index(:,5) = [1, 4, 5];
index(:,4) = [2, 3, 4];
index(:,3) = [2, 3, 5];
index(:,2) = [2, 4, 5];
index(:,1) = [3, 4, 5];

for i = 1:3
	for j = 1:sim_options.num_det2x2
        [Mult_DetM_3x3_int(j,i), Mult_DetM_3x3_overflow(j,i), s.Mult_DetM_3x3_array{j,i}, s.Mult_DetM_3x3_array_mult_total_width{j,i}] = mult(cast(a_int(3,index(i,j)),sim_options.type_3x3_det), cast(DetM_2x2_int(j),sim_options.type_3x3_det), sim_options.type_3x3_det, sim_options.width_hilbert);

        % Проверка переполнения умножителя
        if (Mult_DetM_3x3_overflow(j,i) == 1)
		    disp('Mult in function det3x3 overflow');
            disp({j,i});
            disp({Mult_DetM_3x3_int(j,i)});
        end
	    % Проверка выходной разрядности умножителя
        if (Mult_DetM_3x3_array_mult_total_width(j,i) > sim_options.width_hilbert-1)
		    disp('Mult total width in function det3x3 higher than 64');
            disp({Mult_DetM_3x3_int(j,i)});
        end

        % Наложение маски на первый умножитель
        if sim_options.enable_mask == true
            c1 = bitmask(Mult_DetM_3x3_int(j,i), sim_options.type_3x3_det, width_mult1(i));
            if c1 ~= Mult_DetM_3x3_int(j,i)
                disp('Bit mask error mult det 3x3');
                % c2 = bitmask(sum(i+1,n), sim_options.int_size, width_sum(i+1));
                disp(width_total_mult1(i));
                disp({c1, Mult_DetM_3x3_int(j,i)});
                disp({sim_options.freq, sim_options.SNR});
            end
            Mult_DetM_3x3_int(j,i) = c1;
        end
    end
end

%% double 3x3
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


%% Сумматоры определителя 3х3
[DetM_3x3_int_sum1(1,1), DetM_3x3_int_sum1_overflow(1,1), DetM_3x3_int_sum1_abs(1,1), DetM_3x3_int_sum1_width_total(1,1)] = ...
    adder(Mult_DetM_3x3_int(1,1), -Mult_DetM_3x3_int(2,1), sim_options.type_3x3_det, sim_options.width_hilbert);

[DetM_3x3_int(1,1), DetM_3x3_int_overflow(1,1), DetM_3x3_int_abs(1,1), DetM_3x3_int_width_total(10,1)] = ...
    adder(Mult_DetM_3x3_int(5,1), DetM_3x3_int_sum1(1,1), sim_options.type_3x3_det, sim_options.width_hilbert);
	
[DetM_3x3_int_sum1(2,1), DetM_3x3_int_sum1_overflow(2,1), DetM_3x3_int_sum1_abs(2,1), DetM_3x3_int_sum1_width_total(2,1)] = ...
    adder(Mult_DetM_3x3_int(1,2), -Mult_DetM_3x3_int(3,1), sim_options.type_3x3_det, sim_options.width_hilbert);
[DetM_3x3_int(2,1), DetM_3x3_int_overflow(2,1), DetM_3x3_int_abs(2,1), DetM_3x3_int_width_total(2,1)] = ...
    adder(Mult_DetM_3x3_int(6,1), DetM_3x3_int_sum1(2,1), sim_options.type_3x3_det, sim_options.width_hilbert);

[DetM_3x3_int_sum1(3,1), DetM_3x3_int_sum1_overflow(3,1), DetM_3x3_int_sum1_abs(3,1), DetM_3x3_int_sum1_width_total(3,1)] = ...
    adder(Mult_DetM_3x3_int(1,3), -Mult_DetM_3x3_int(4,1), sim_options.type_3x3_det, sim_options.width_hilbert);
[DetM_3x3_int(3,1), DetM_3x3_int_overflow(3,1), DetM_3x3_int_abs(3,1), DetM_3x3_int_width_total(3,1)] = ...
    adder(Mult_DetM_3x3_int(7,1), DetM_3x3_int_sum1(3,1), sim_options.type_3x3_det, sim_options.width_hilbert);

[DetM_3x3_int_sum1(4,1), DetM_3x3_int_sum1_overflow(4,1), DetM_3x3_int_sum1_abs(4,1), DetM_3x3_int_sum1_width_total(4,1)] = ...
    adder(Mult_DetM_3x3_int(2,2), -Mult_DetM_3x3_int(3,2), sim_options.type_3x3_det, sim_options.width_hilbert);
[DetM_3x3_int(4,1), DetM_3x3_int_overflow(4,1), DetM_3x3_int_abs(4,1), DetM_3x3_int_width_total(4,1)] = ...
    adder(Mult_DetM_3x3_int(8,1), DetM_3x3_int_sum1(4,1), sim_options.type_3x3_det, sim_options.width_hilbert);

[DetM_3x3_int_sum1(5,1), DetM_3x3_int_sum1_overflow(5,1), DetM_3x3_int_sum1_abs(5,1), DetM_3x3_int_sum1_width_total(5,1)] = ...
    adder(Mult_DetM_3x3_int(2,3), -Mult_DetM_3x3_int(4,2), sim_options.type_3x3_det, sim_options.width_hilbert);
[DetM_3x3_int(5,1), DetM_3x3_int_overflow(5,1), DetM_3x3_int_abs(5,1), DetM_3x3_int_width_total(5,1)] = ...
    adder(Mult_DetM_3x3_int(9,1), DetM_3x3_int_sum1(5,1), sim_options.type_3x3_det, sim_options.width_hilbert);

[DetM_3x3_int_sum1(6,1), DetM_3x3_int_sum1_overflow(6,1), DetM_3x3_int_sum1_abs(6,1), DetM_3x3_int_sum1_width_total(6,1)] = ...
    adder(Mult_DetM_3x3_int(3,3), -Mult_DetM_3x3_int(4,3), sim_options.type_3x3_det, sim_options.width_hilbert);
[DetM_3x3_int(6,1), DetM_3x3_int_overflow(6,1), DetM_3x3_int_abs(6,1), DetM_3x3_int_width_total(6,1)] = ...
    adder(Mult_DetM_3x3_int(10,1), DetM_3x3_int_sum1(6,1), sim_options.type_3x3_det, sim_options.width_hilbert);

[DetM_3x3_int_sum1(7,1), DetM_3x3_int_sum1_overflow(7,1), DetM_3x3_int_sum1_abs(7,1), DetM_3x3_int_sum1_width_total(7,1)] = ...
    adder(Mult_DetM_3x3_int(5,2), -Mult_DetM_3x3_int(6,2), sim_options.type_3x3_det, sim_options.width_hilbert);
[DetM_3x3_int(7,1), DetM_3x3_int_overflow(7,1), DetM_3x3_int_abs(7,1), DetM_3x3_int_width_total(7,1)] = ...
    adder(Mult_DetM_3x3_int(8,2), DetM_3x3_int_sum1(7,1), sim_options.type_3x3_det, sim_options.width_hilbert);

[DetM_3x3_int_sum1(8,1), DetM_3x3_int_sum1_overflow(8,1), DetM_3x3_int_sum1_abs(8,1), DetM_3x3_int_sum1_width_total(8,1)] = ...
    adder(Mult_DetM_3x3_int(5,3), -Mult_DetM_3x3_int(7,2), sim_options.type_3x3_det, sim_options.width_hilbert);
[DetM_3x3_int(8,1), DetM_3x3_int_overflow(8,1), DetM_3x3_int_abs(8,1), DetM_3x3_int_width_total(3,1)] = ...
    adder(Mult_DetM_3x3_int(9,2), DetM_3x3_int_sum1(8,1), sim_options.type_3x3_det, sim_options.width_hilbert);

[DetM_3x3_int_sum1(9,1), DetM_3x3_int_sum1_overflow(9,1), DetM_3x3_int_sum1_abs(9,1), DetM_3x3_int_sum1_width_total(9,1)] = ...
    adder(Mult_DetM_3x3_int(6,3), -Mult_DetM_3x3_int(7,3), sim_options.type_3x3_det, sim_options.width_hilbert);
[DetM_3x3_int(9,1), DetM_3x3_int_overflow(9,1), DetM_3x3_int_abs(9,1), DetM_3x3_int_width_total(9,1)] = ...
    adder(Mult_DetM_3x3_int(10,2), DetM_3x3_int_sum1(9,1), sim_options.type_3x3_det, sim_options.width_hilbert);

[DetM_3x3_int_sum1(10,1), DetM_3x3_int_sum1_overflow(10,1), DetM_3x3_int_sum1_abs(10,1), DetM_3x3_int_sum1_width_total(10,1)] = ...
    adder(Mult_DetM_3x3_int(8,3), -Mult_DetM_3x3_int(9,3), sim_options.type_3x3_det, sim_options.width_hilbert);
[DetM_3x3_int(10,1), DetM_3x3_int_overflow(10,1), DetM_3x3_int_abs(10,1), DetM_3x3_int_width_total(10,1)] = ...
    adder(Mult_DetM_3x3_int(10,3), DetM_3x3_int_sum1(10,1), sim_options.type_3x3_det, sim_options.width_hilbert);


for i = 1:sim_options.num_det2x2
    % Проверка переполнения пресумматора
    if (DetM_3x3_int_sum1_overflow(i,1) == 1)
		disp('Sum det3x3 overflow');
        disp({DetM_3x3_int_sum1(i,1)});
        disp({i});
    end
	% Проверка выходной разрядности пресумматора
    if (DetM_3x3_int_sum1_width_total(1,1) > sim_options.width_hilbert-1)
		disp('Sum total width in function det3x3 higher than 64');
        disp({DetM_3x3_int_sum1(i,1)});
        disp({i});
    end
    % Наложение маски на пресумматор
    if sim_options.enable_mask == true
        c1 = bitmask(DetM_3x3_int_sum1(i,1), sim_options.type_3x3_det, width_mult1(i));
        if c1 ~= DetM_3x3_int_sum1(i,i)
            disp('Bit mask error pre_sum det 3x3');
            % c2 = bitmask(sum(i+1,n), sim_options.int_size, width_sum(i+1));
            disp(DetM_3x3_int_sum1_width_total(i));
            disp({c1, DetM_3x3_int_sum1(i,1)});
            disp({sim_options.freq, sim_options.SNR});
        end
        DetM_3x3_int_sum1(i,1) = c1;
    end

    %% Проверка переполнения сумматора
    if (DetM_3x3_int_overflow(i,1) == 1)
		disp('Sum det3x3 overflow');
        disp({DetM_3x3_int(i,1)});
        disp({i});
    end
	% Проверка выходной разрядности сумматора
    if (DetM_3x3_int_width_total(i,1) > sim_options.width_hilbert-1)
		disp('Sum total width in function det3x3 higher than 64');
        disp({DetM_3x3_int(i,1)});
        disp({i});
    end
    % Наложение маски на сумматор
    if sim_options.enable_mask == true
        c1 = bitmask(DetM_3x3_int(i,1), sim_options.type_3x3_det, width_mult1(i));
        if c1 ~= DetM_3x3_int(i,i)
            disp('Bit mask error pre_sum det 3x3');
            % c2 = bitmask(sum(i+1,n), sim_options.int_size, width_sum(i+1));
            disp(DetM_3x3_int_width_total(i));
            disp({c1, DetM_3x3_int(i,1)});
            disp({sim_options.freq, sim_options.SNR});
        end
        DetM_3x3_int(i,1) = c1;
    end

    s.DetM_3x3_int_pre_sum_array{i} = DetM_3x3_int_sum1_abs(i,1);
    s.DetM_3x3_int_pre_sum_width_total{i} = DetM_3x3_int_sum1_width_total(i,1);
    s.DetM_3x3_int_sum_array{i} = DetM_3x3_int_abs(i,1);
    s.DetM_3x3_int_sum_width_total{i} = DetM_3x3_int_width_total(i,1);
end

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

%% Умножители определителя 4х4
index1(1,:) = [2, 3, 4, 5];
index1(2,:) = [1, 3, 4, 5];
index1(3,:) = [1, 2, 4, 5];
index1(4,:) = [1, 2, 3, 5];
index1(5,:) = [1, 2, 3, 4];

index_DetM_3x3(1,:) = [10, 9, 8, 7];
index_DetM_3x3(2,:) = [10, 6, 5, 4];
index_DetM_3x3(3,:) = [9, 6, 3, 2];
index_DetM_3x3(4,:) = [8, 5, 3, 1]; 
index_DetM_3x3(5,:) = [7, 4, 2, 1];

kk = 0;
for j = 1:5
    for i = 1:4
        [var_mult_det4x4_int(j,i), var_mult_det4x4_int_overflow(j,i), var_mult_det4x4_abs(j,i), var_mult_det4x4_width_total(j,i)] = mult(cast(a_int(2,index1(j,i)), sim_options.type_4x4_det), ...
            cast(DetM_3x3_int(index_DetM_3x3(j,i)), sim_options.type_4x4_det), sim_options.type_4x4_det, sim_options.width_hilbert);

        % Проверка переполнения умножителя
        if (var_mult_det4x4_int_overflow(j,i) == 1)
		    disp('Mult det4x4 overflow');
            disp({var_mult_det4x4_int(j,i)});
            disp({j,i});
        end
	    % Проверка выходной разрядности умножителя
        if (var_mult_det4x4_width_total(j,i) > sim_options.width_hilbert-1)
		    disp('Mult total width in function det4x4 higher than 64');
            disp({var_mult_det4x4_int(i,1)});
            disp({i});
        end
        % Наложение маски на умножитель
        if sim_options.enable_mask == true
            c1 = bitmask(var_mult_det4x4_int(j,i), sim_options.type_3x3_det, width_mult1(i));
            if c1 ~= var_mult_det4x4_int(j,i)
                disp('Bit mask error pre_sum det 3x3');
                % c2 = bitmask(sum(i+1,n), sim_options.int_size, width_sum(i+1));
                disp(var_mult_det4x4_width_total(j,i));
                disp({c1, var_mult_det4x4_int(j,i)});
                disp({sim_options.freq, sim_options.SNR});
            end
            var_mult_det4x4_int(j,i) = c1;
        end        
    end
    DetM_4x4_int_mult_array(kk+1:kk+4) = var_mult_det4x4_abs(j,:);
    DetM_4x4_int_mult_width_total(kk+1:kk+4) = var_mult_det4x4_width_total(j,:);
    kk = kk + 4;
end

%% Сумматоры определителя 4х4
kk = 0;
for i = 1:5
    for j = 1:2
        if (mod(j,2) == 1)
            [DetM_4x4_int_sum(i,j), DetM_4x4_int_sum1_overflow(i,j), DetM_4x4_int_sum1_abs(i,j), DetM_4x4_int_sum1_width_total(i,j)] = adder(var_mult_det4x4_int(i,j), ...
                -var_mult_det4x4_int(i,j+1),  sim_options.type_4x4_det, sim_options.width_hilbert);
        
            % Проверка переполнения сумматора
            if (DetM_4x4_int_sum1_overflow(i,j) == 1)
		        disp('Sum1 det4x4 overflow');
                disp({DetM_4x4_int_sum(i,j)});
                disp({i,j});
            end
	        % Проверка выходной разрядности сумматора
            if (DetM_4x4_int_sum1_width_total(i,j) > sim_options.width_hilbert-1)
		        disp('Sum1 total width in function det4x4 higher than 64');
                disp({DetM_4x4_int_sum(i,1)});
                disp({i,j});
            end

            % Наложение маски на сумматор
            if sim_options.enable_mask == true
                c1 = bitmask(DetM_4x4_int_sum(i,j), sim_options.type_4x4_det, width_mult1(i));
                if c1 ~= DetM_4x4_int_sum(i,j)
                    disp('Bit mask error pre_sum det 4x4');
                    % c2 = bitmask(sum(i+1,n), sim_options.int_size, width_sum(i+1));
                    % disp(DetM_4x4_int_sum(i,j));
                    disp({c1, DetM_4x4_int_sum(i,j)});
                    disp({sim_options.freq, sim_options.SNR});
                end
                DetM_4x4_int_sum(i,j) = c1;
            end
        else
            [DetM_4x4_int_sum(i,j), DetM_4x4_int_sum1_overflow(i,j), DetM_4x4_int_sum1_abs(i,j), DetM_4x4_int_sum1_width_total(i,j)] = adder(var_mult_det4x4_int(i,j+1), ...
                -var_mult_det4x4_int(i,j+2),  sim_options.type_4x4_det, sim_options.width_hilbert);
         
            % Проверка переполнения сумматора
            if (DetM_4x4_int_sum1_overflow(i,j) == 1)
		        disp('Sum2 det4x4 overflow');
                disp({DetM_4x4_int_sum(i,j)});
                disp({i,j});
            end
	        % Проверка выходной разрядности сумматора
            if (DetM_4x4_int_sum1_width_total(i,j) > sim_options.width_hilbert-1)
		        disp('Sum2 total width in function det4x4 higher than 64');
                disp({DetM_4x4_int_sum(i,j)});
                disp({i,j});
            end

            % Наложение маски на сумматор
            if sim_options.enable_mask == true
                c1 = bitmask(DetM_4x4_int_sum(i,j), sim_options.type_4x4_det, width_mult1(i));
                if c1 ~= DetM_4x4_int_sum(i,j)
                    disp('Bit mask error pre_sum det 4x4');
                    % c2 = bitmask(sum(i+1,n), sim_options.int_size, width_sum(i+1));
                    % disp(DetM_4x4_int_sum(i,j));
                    disp({c1, DetM_4x4_int_sum(i,j)});
                    disp({sim_options.freq, sim_options.SNR});
                end
                DetM_4x4_int_sum(i,j) = c1;
            end
        end
        DetM_4x4_int_pre_sum_array(kk+1) = DetM_4x4_int_sum1_abs(i,j);
        DetM_4x4_int_pre_sum_width_total(kk+1) = DetM_4x4_int_sum1_width_total(i,j);
        kk = kk + 1;
    end
    [DetM_4x4_int(i), DetM_4x4_int_overflow(i), DetM_4x4_int_abs(i), DetM_4x4_int_width_total(i)] = adder(DetM_4x4_int_sum(i,1), ...
        DetM_4x4_int_sum(i,2), sim_options.type_4x4_det, sim_options.width_hilbert);

    % Проверка переполнения сумматора
    if (DetM_4x4_int_overflow(i) == 1)
	    disp('Sum det4x4 overflow');
        disp({DetM_4x4_int(i)});
        disp({i});
    end
	% Проверка выходной разрядности сумматора
    if (DetM_4x4_int_width_total(i) > sim_options.width_hilbert-1)
	    disp('Sum total width in function det4x4 higher than 64');
        disp({DetM_4x4_int(i)});
        disp({i});
    end

    % Наложение маски на сумматор
    if sim_options.enable_mask == true
        c1 = bitmask(DetM_4x4_int(i), sim_options.type_4x4_det, width_mult1(i));
        if c1 ~= DetM_4x4_int(i)
            disp('Bit mask error pre_sum det 3x3');
            % c2 = bitmask(sum(i+1,n), sim_options.int_size, width_sum(i+1));
            % disp(DetM_4x4_int_sum(i,j));
            disp({c1, DetM_4x4_int(i)});
            disp({sim_options.freq, sim_options.SNR});
        end
        DetM_4x4_int(i) = c1;
    end

    DetM_4x4_int_sum_array(i) = DetM_4x4_int_abs(i);
    DetM_4x4_int_sum_array_width_total(i) = DetM_4x4_int_width_total(i);
    
end

                                                                            %% 5x5
	% double                                                                    
	var1 = a(1,1) * DetM_4x4_n_1;
	var2 = a(1,2) * DetM_4x4_n_2;
	var3 = a(1,3) * DetM_4x4_n_3;
	var4 = a(1,4) * DetM_4x4_n_4;
	var5 = a(1,5) * DetM_4x4_n_5;
	DetM_5x5 = var1 - var2 + var3 - var4 + var5;

	%% Умножители матрицы 5х5
	for i = 1:5
		[var_int_det5x5(i,:), var_int_det5x5_overflow(i,:), var_int_det5x5_abs(i,:), var_int_det5x5_mult_width_total(i,:)] = mult(cast(a_int(1,i),sim_options.type_5x5_det), ...
			DetM_4x4_int(i), sim_options.type_5x5_det, sim_options.width_hilbert);
	
		% Проверка переполнения сумматора
		if (var_int_det5x5_overflow(i) == 1)
			disp('Mult det5x5 overflow');
			disp({var_int_det5x5(i)});
			disp({i});
		end
		% Проверка выходной разрядности сумматора
		if (var_int_det5x5_mult_width_total(i) > sim_options.width_hilbert-1)
			disp('Mult total width in function det5x5 higher than 64');
			disp({var_int_det5x5(i)});
			disp({i});
		end

		% Наложение маски на сумматор
		if sim_options.enable_mask == true
			c1 = bitmask(var_int_det5x5(i), sim_options.type_5x5_det, width_mult1(i));
			if c1 ~= var_int_det5x5(i)
				disp('Bit mask error mult det 5x5');
				% c2 = bitmask(sum(i+1,n), sim_options.int_size, width_sum(i+1));
				% disp(DetM_4x4_int_sum(i,j));
				disp({c1, var_int_det5x5(i)});
				disp({sim_options.freq, sim_options.SNR});
			end
			var_int_det5x5(i) = c1;
		end
	
		DetM_5x5_int_mult_array(i) = var_int_det5x5_abs(i);
		DetM_5x5_int_mult_array_width_total(i) = var_int_det5x5_mult_width_total(i);
	
	end

	%% Сумматоры1 определителя 5х5
	for i = 1:2
		if (mod(i,2) == 1)
			[DetM_5x5_int_sum1(i), DetM_5x5_int_sum1_overflow(i), DetM_5x5_int_sum1_abs(i), DetM_5x5_int_sum1_width_total(i)] = adder(var_int_det5x5(i), ...
				-var_int_det5x5(i+1), sim_options.type_5x5_det, sim_options.width_hilbert);
			% Проверка переполнения сумматора
			if (DetM_5x5_int_sum1_overflow(i) == 1)
				disp('Sum1 det5x5 overflow');
				disp({DetM_5x5_int_sum1_overflow(i)});
				disp({i});
			end
			% Проверка выходной разрядности сумматора
			if (DetM_5x5_int_sum1_width_total(i) > sim_options.width_hilbert-1)
				disp('Sum1 total width in function det5x5 higher than 64');
				disp({DetM_5x5_int_sum1_width_total(i)});
				disp({i});
			end
			% Наложение маски на сумматор
			if sim_options.enable_mask == true
				c1 = bitmask(DetM_5x5_int_sum1(i), sim_options.type_5x5_det, width_mult1(i));
				if c1 ~= DetM_5x5_int_sum1(i)
					disp('Bit mask error mult det 5x5');
					% c2 = bitmask(sum(i+1,n), sim_options.int_size, width_sum(i+1));
					% disp(DetM_4x4_int_sum(i,j));
					disp({c1, DetM_5x5_int_sum1(i)});
					disp({sim_options.freq, sim_options.SNR});
				end
				DetM_5x5_int_sum1(i) = c1;
			end
		else
			[DetM_5x5_int_sum1(i), DetM_5x5_int_sum1_overflow(i), DetM_5x5_int_sum1_abs(i), DetM_5x5_int_sum1_width_total(i)] = adder(var_int_det5x5(i+1), ...
				-var_int_det5x5(i+2), sim_options.type_5x5_det, sim_options.width_hilbert);
			% Проверка переполнения сумматора
			if (DetM_5x5_int_sum1_overflow(i) == 1)
				disp('Sum1 det5x5 overflow');
				disp({DetM_5x5_int_sum1(i)});
				disp({i});
			end
			% Проверка выходной разрядности сумматора
			if (DetM_5x5_int_sum1_width_total(i) > sim_options.width_hilbert-1)
				disp('Sum1 total width in function det5x5 higher than 64');
				disp({DetM_5x5_int_sum1(i)});
				disp({i});
			end
			% Наложение маски на сумматор
			if sim_options.enable_mask == true
				c1 = bitmask(DetM_5x5_int_sum1(i), sim_options.type_5x5_det, width_mult1(i));
				if c1 ~= DetM_5x5_int_sum1(i)
					disp('Bit mask error mult det 5x5');
					% c2 = bitmask(sum(i+1,n), sim_options.int_size, width_sum(i+1));
					% disp(DetM_4x4_int_sum(i,j));
					disp({c1, DetM_5x5_int_sum1(i)});
					disp({sim_options.freq, sim_options.SNR});
				end
				DetM_5x5_int_sum1(i) = c1;
			end
		end
		DetM_5x5_int_pre_sum1_array(i) = DetM_5x5_int_sum1_abs(i);
		DetM_5x5_int_pre_sum1_array_width_total(i) = DetM_5x5_int_sum1_width_total(i);
	end

	%% Сумматоры2 определителя 5х5
	[DetM_5x5_int_sum3, DetM_5x5_int_sum3_overflow, DetM_5x5_int_sum3_abs, DetM_5x5_int_sum3_width_total] = adder(DetM_5x5_int_sum1(1), ...
		DetM_5x5_int_sum1(2), 	sim_options.type_5x5_det, sim_options.width_hilbert);
	
	% Проверка переполнения сумматора
	if (DetM_5x5_int_sum3_overflow == 1)
		disp('Sum3 det5x5 overflow');
		disp({DetM_5x5_int_sum3});
	end
	% Проверка выходной разрядности сумматора
	if (DetM_5x5_int_sum3_width_total > sim_options.width_hilbert-1)
		disp('Sum3 total width in function det5x5 higher than 64');
		disp({DetM_5x5_int_sum3});
	end

	% Наложение маски на сумматор
	if sim_options.enable_mask == true
		c1 = bitmask(DetM_5x5_int_sum3, sim_options.type_5x5_det, width_mult1(i));
		if c1 ~= DetM_5x5_int_sum3
			disp('Bit mask error mult det 5x5');
			% c2 = bitmask(sum(i+1,n), sim_options.int_size, width_sum(i+1));
			% disp(DetM_4x4_int_sum(i,j));
			disp({c1, DetM_5x5_int_sum3});
			disp({sim_options.freq, sim_options.SNR});
		end
		DetM_5x5_int_sum3 = c1;
	end

	%% Сумматоры3 определителя 5х5
	[DetM_5x5_int, DetM_5x5_int_overflow, DetM_5x5_int_abs, DetM_5x5_int_width_total] = adder(DetM_5x5_int_sum3,  var_int_det5x5(5), ...
		sim_options.type_5x5_det, sim_options.width_hilbert);

	% Проверка переполнения сумматора
	if (DetM_5x5_int_overflow == 1)
		disp('Sum3 det5x5 overflow');
		disp({DetM_5x5_int});
	end
	% Проверка выходной разрядности сумматора
	if (DetM_5x5_int_width_total > sim_options.width_hilbert-1)
		disp('Sum3 total width in function det5x5 higher than 64');
		disp({DetM_5x5_int});
	end
	% Наложение маски на сумматор
	if sim_options.enable_mask == true
		c1 = bitmask(DetM_5x5_int, sim_options.type_5x5_det, width_mult1(i));
		if c1 ~= DetM_5x5_int
			disp('Bit mask error mult det 5x5');
			% c2 = bitmask(sum(i+1,n), sim_options.int_size, width_sum(i+1));
			% disp(DetM_4x4_int_sum(i,j));
			disp({c1, DetM_5x5_int});
			disp({sim_options.freq, sim_options.SNR});
		end
		DetM_5x5_int = c1;
	end





















        
	%%
	for n = 1:sim_options.num_det2x2
		if Det_3x3_LU_matlab(n) < 0
			Det_3x3_LU_matlab(n) = Det_3x3_LU_matlab(n) * -1;
		end
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
	
	for n = 1:5
		if Det_4x4_LU_matlab(n) < 0
			Det_4x4_LU_matlab(n) = Det_4x4_LU_matlab(n) * -1;
		end
		if DetM_4x4_array(n) < 0
			DetM_4x4_array(n) = DetM_4x4_array(n) * -1;
		end
	end
end