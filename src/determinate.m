
function [DetM_5x5, DetM_5x5_int, DetM_2x2_abs, Det_2x2_LU_matlab, DetM_3x3_array, Det_3x3_LU_matlab, DetM_4x4_array, Det_4x4_LU_matlab, s ...
    ] = determinate(data_in, data_in_int, det_max_width, sim_options)

s = struct;

DetM_2x2 = zeros(sim_options.num_det2x2,1);
DetM_2x2_int = cast(zeros(sim_options.num_det2x2,1), sim_options.type_2x2_det);

Det_2x2_LU_matlab = zeros(sim_options.num_det2x2,1);

s.Det2x2_mult_abs = cast(zeros(sim_options.num_det2x2,1), sim_options.type_2x2_det);
s.Det2x2_sum_abs = cast(zeros(sim_options.num_det2x2,1), sim_options.type_2x2_det);

%% Определитель 3х3
Mult_DetM_3x3_int = cast(zeros(sim_options.num_det2x2,3), sim_options.type_3x3_det);

%%
Det_3x3_LU_matlab 							= zeros(sim_options.num_det2x2,1);
Mult_DetM_3x3_array                         = cast(zeros(sim_options.num_det2x2,3), sim_options.type_3x3_det);
Mult_DetM_3x3_array_mult_total_width 		= cast(zeros(sim_options.num_det2x2,3), sim_options.type_3x3_det);
s.Mult_DetM_3x3_array 						= cast(zeros(sim_options.num_det2x2*3,1), sim_options.type_3x3_det);
s.Mult_DetM_3x3_array_mult_total_width 	    = cast(zeros(sim_options.num_det2x2*3,1), sim_options.type_3x3_det);
s.DetM_3x3_int_pre_sum_array 				= cast(zeros(sim_options.num_det2x2,1), sim_options.type_3x3_det);
s.DetM_3x3_int_pre_sum_width_total 		    = cast(zeros(sim_options.num_det2x2,1), sim_options.type_3x3_det);
s.DetM_3x3_int_sum_array 					= cast(zeros(sim_options.num_det2x2,1), sim_options.type_3x3_det);
s.DetM_3x3_int_sum_width_total 			    = cast(zeros(sim_options.num_det2x2,1), sim_options.type_3x3_det);

s.DetM_4x4_int_mult_array	 				= cast(zeros(sim_options.num_det2x2*2,1), sim_options.type_4x4_det);
s.DetM_4x4_int_mult_width_total 			= cast(zeros(sim_options.num_det2x2*2,1), sim_options.type_4x4_det);
s.DetM_4x4_int_pre_sum_array 				= cast(zeros(sim_options.num_det2x2,1), sim_options.type_4x4_det);
s.DetM_4x4_int_pre_sum_width_total 		    = cast(zeros(sim_options.num_det2x2,1), sim_options.type_4x4_det);
s.DetM_4x4_int_sum_array 					= cast(zeros(5,1), sim_options.type_4x4_det);
s.DetM_4x4_int_sum_array_width_total 		= cast(zeros(5,1), sim_options.type_4x4_det);

s.DetM_5x5_int_mult_array 					= cast(zeros(5,1), sim_options.type_5x5_det);
s.DetM_5x5_int_mult_array_width_total 		= cast(zeros(5,1), sim_options.type_5x5_det);
s.DetM_5x5_int_pre_sum1_array 				= cast(zeros(2,1), sim_options.type_5x5_det);
s.DetM_5x5_int_pre_sum1_array_width_total 	= cast(zeros(2,1), sim_options.type_5x5_det);

Mult_DetM_3x3_array_int = cast(zeros(sim_options.num_det2x2*3,1), sim_options.type_3x3_det);
                                                                            %% 2x2
%% Находим матрицы 2x2
e = 0;
t = 4;
for n = 1:4
    for i = 1:t
        for j = 1:1
            s.a{i+e}(:,j) = data_in(4:end,n);
            s.a_int{i+e}(:,j) = data_in_int(4:end,n);
        end
        s.a{i+e}(:,2) = data_in(4:end,n+i);
        s.a_int{i+e}(:,2) = data_in_int(4:end,n+i);
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

    %% Определитель 2х2 double
    DetM_2x2(i) = s.a{i}(1,1) * s.a{i}(2,2) - s.a{i}(2,1) * s.a{i}(1,2);

    DetM_2x2_abs(i) = DetM_2x2(i);
    if DetM_2x2_abs(i) < 0
        DetM_2x2_abs(i) = DetM_2x2_abs(i) * -1;
    end

    %% Умножители определителя 2х2
    [mult_int(i+(i-1)), mult_overflow(i+(i-1)), s.Det2x2_mult_abs(i+(i-1)), width_total_mult(i+(i-1))] = mult(cast(s.a_int{i}(1,1),sim_options.type_2x2_det), ...
        cast(s.a_int{i}(2,2),sim_options.type_2x2_det), sim_options.type_2x2_det, sim_options.width_fractional);

    %% Проверка переполнения умножителя
    if (mult_overflow(i+(i-1)) == 1)
	    disp(['Mult in function det2x2 overflow in ', num2str(i+(i-1))]);
        disp(sim_options.width_hilbert);
    end
	% Проверка выходной разрядности умножителя
    if (width_total_mult(i+(i-1)) > sim_options.width_hilbert)
	    disp(['Mult total width in function det2x2 higher than 64 ', num2str(i+(i-1))]);
        disp(sim_options.width_hilbert);
    end

    % Наложение маски на первый умножитель
    if sim_options.enable_mask == true
        c1 = bitmask(mult_int(i+(i-1)), sim_options.type_2x2_det, det_max_width(i+(i-1),1));
        if c1 ~= mult_int(i+(i-1))
            disp(['Bit mask error mult det 2x2 in ', num2str(i+(i-1))]);
            % c2 = bitmask(sum(i+1,n), sim_options.int_size, width_sum(i+1));
            disp(width_total_mult(i+(i-1)));
            disp([c1, mult_int(i+(i-1))]);
            disp([sim_options.freq, sim_options.SNR]);
        end
        mult_int(i+(i-1)) = c1;
    end

    [mult_int(i+i), mult_overflow(i+i), s.Det2x2_mult_abs(i+i), width_total_mult(i+i)] = mult(cast(s.a_int{i}(2,1),sim_options.type_2x2_det), ...
        cast(s.a_int{i}(1,2),sim_options.type_2x2_det), sim_options.type_2x2_det, sim_options.width_fractional);

    %% Проверка переполнения умножителя
    if (mult_overflow(i+i) == 1)
		disp(['Mult in function det2x2 overflow in ', num2str(i+i)]);
        disp(sim_options.width_hilbert);
    end
	% Проверка выходной разрядности умножителя
    if (width_total_mult(i+i) > sim_options.width_hilbert)
		disp(['Mult total width in function det2x2 higher than 64 ', num2str(i+i)]);
        disp(sim_options.width_hilbert);
    end

    % Наложение маски на второй умножитель
    if sim_options.enable_mask == true
        c1 = bitmask(mult_int(i+i), sim_options.type_2x2_det, det_max_width(i+i,1));
        if c1 ~= mult_int(i+i)
            disp(['Bit mask error mult 2 det 2x2 in ', num2str(i+i)]);
            disp(width_total_mult(i+i));
            disp([c1, mult_int(i+i)]);
            disp([sim_options.freq, sim_options.SNR]);
        end
        mult_int(i+i) = c1;
    end

    %% Сумматоры определителя 2х2
	[DetM_2x2_int(i), sum_overflow(i), s.Det2x2_sum_abs(i), width_total_sum(i)] = adder(cast(mult_int(i+(i-1)),sim_options.type_2x2_det), ...
        cast(-mult_int(i+i),sim_options.type_2x2_det), sim_options.type_2x2_det, sim_options.width_fractional);

    % Проверка переполнения сумматора
    if (sum_overflow(i) == 1)
		disp(['Adder in function det2x2 overflow in ', num2str(i)]);
        disp(sim_options.width_hilbert);
        disp([s.a_int{i}]);
    end
	% Проверка выходной разрядности сумматора
    if (width_total_sum(i) > sim_options.width_hilbert)
		disp(['Sum total width in function det2x2 higher than 64 in ', num2str(i)]);
        disp(sim_options.width_hilbert);
        disp([width_total_sum(i)]);
    end
    % Наложение маски на сумматор
    if sim_options.enable_mask == true
        c1 = bitmask(DetM_2x2_int(i), sim_options.type_2x2_det, det_max_width(i,2));
        if c1 ~= DetM_2x2_int(i)
            disp(['Bit mask error sum det 2x2 in ', num2str(i)]);
            disp([width_total_sum(i), det_max_width(i,2)]);
            disp([c1, DetM_2x2_int(i)]);
            disp([sim_options.freq, sim_options.SNR]);
        end
        DetM_2x2_int(i) = c1;
    end
end
                                                                            %% 3x3
% double
Mult_DetM_2x2_n_10_a31 = data_in(3,1) * DetM_2x2(10);
Mult_DetM_2x2_n_10_a32 = data_in(3,2) * DetM_2x2(10);
Mult_DetM_2x2_n_10_a33 = data_in(3,3) * DetM_2x2(10);

Mult_DetM_2x2_n_9_a31 = data_in(3,1) * DetM_2x2(9);
Mult_DetM_2x2_n_9_a32 = data_in(3,2) * DetM_2x2(9);
Mult_DetM_2x2_n_9_a34 = data_in(3,4) * DetM_2x2(9);

Mult_DetM_2x2_n_8_a31 = data_in(3,1) * DetM_2x2(8);
Mult_DetM_2x2_n_8_a32 = data_in(3,2) * DetM_2x2(8);
Mult_DetM_2x2_n_8_a35 = data_in(3,5) * DetM_2x2(8);

Mult_DetM_2x2_n_7_a31 = data_in(3,1) * DetM_2x2(7);
Mult_DetM_2x2_n_7_a33 = data_in(3,3) * DetM_2x2(7);
Mult_DetM_2x2_n_7_a34 = data_in(3,4) * DetM_2x2(7);

Mult_DetM_2x2_n_6_a31 = data_in(3,1) * DetM_2x2(6);
Mult_DetM_2x2_n_6_a33 = data_in(3,3) * DetM_2x2(6);
Mult_DetM_2x2_n_6_a35 = data_in(3,5) * DetM_2x2(6);

Mult_DetM_2x2_n_5_a31 = data_in(3,1) * DetM_2x2(5);
Mult_DetM_2x2_n_5_a34 = data_in(3,4) * DetM_2x2(5);
Mult_DetM_2x2_n_5_a35 = data_in(3,5) * DetM_2x2(5);

Mult_DetM_2x2_n_4_a32 = data_in(3,2) * DetM_2x2(4);
Mult_DetM_2x2_n_4_a33 = data_in(3,3) * DetM_2x2(4);
Mult_DetM_2x2_n_4_a34 = data_in(3,4) * DetM_2x2(4);

Mult_DetM_2x2_n_3_a32 = data_in(3,2) * DetM_2x2(3);
Mult_DetM_2x2_n_3_a33 = data_in(3,3) * DetM_2x2(3);
Mult_DetM_2x2_n_3_a35 = data_in(3,5) * DetM_2x2(3);

Mult_DetM_2x2_n_2_a32 = data_in(3,2) * DetM_2x2(2);
Mult_DetM_2x2_n_2_a34 = data_in(3,4) * DetM_2x2(2);
Mult_DetM_2x2_n_2_a35 = data_in(3,5) * DetM_2x2(2);

Mult_DetM_2x2_n_1_a33 = data_in(3,3) * DetM_2x2(1);
Mult_DetM_2x2_n_1_a34 = data_in(3,4) * DetM_2x2(1);
Mult_DetM_2x2_n_1_a35 = data_in(3,5) * DetM_2x2(1);

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
        [Mult_DetM_3x3_int(j,i), Mult_DetM_3x3_overflow(j,i), Mult_DetM_3x3_array(j,i), Mult_DetM_3x3_array_mult_total_width(j,i)] = mult(cast(data_in_int(3,index(i,j)),sim_options.type_3x3_det), cast(DetM_2x2_int(j),sim_options.type_3x3_det), sim_options.type_3x3_det, sim_options.width_hilbert);

        % Проверка переполнения умножителя
        if (Mult_DetM_3x3_overflow(j,i) == 1)
		    disp(['Mult in function det3x3 overflow in ', num2str(j), num2str(i)]);
            disp({Mult_DetM_3x3_int(j,i)});
        end
	    % Проверка выходной разрядности умножителя
        if (Mult_DetM_3x3_array_mult_total_width(j,i) > sim_options.width_hilbert-1)
		    disp(['Mult total width in function det3x3 higher than 64 in ', num2str(j), num2str(i)]);
            disp([Mult_DetM_3x3_int(j,i)]);
        end
    end
    s.Mult_DetM_3x3_array(i:3:end) = Mult_DetM_3x3_array(:,i);
    s.Mult_DetM_3x3_array_mult_total_width(i:3:end) = Mult_DetM_3x3_array_mult_total_width(:,i);

    Mult_DetM_3x3_array_int(i:3:end) = Mult_DetM_3x3_int(:,i);
end

for j = 1:sim_options.num_det2x2*3
    % Наложение маски на первый умножитель определителя 3х3
    if sim_options.enable_mask == true
        c1 = bitmask(Mult_DetM_3x3_array_int(j), sim_options.type_3x3_det, det_max_width(j,3));
        if c1 ~= Mult_DetM_3x3_array_int(j)
            disp(['Bit mask error mult det 3x3 in ', num2str(j)]);
            disp([c1, Mult_DetM_3x3_array_int(j)]);
            disp([sim_options.freq, sim_options.SNR]);
        end
        Mult_DetM_3x3_array_int(j) = c1;
    end
end


%% double 3x3
DetM_3x3_n_11 = Mult_DetM_2x2_n_10_a33 - Mult_DetM_2x2_n_9_a34 + Mult_DetM_2x2_n_8_a35; 
Det_3x3_LU_matlab(1) = (det(data_in(3:end,3:end))); % 3 4 5

DetM_3x3_n_12 = Mult_DetM_2x2_n_10_a32 - Mult_DetM_2x2_n_7_a34 + Mult_DetM_2x2_n_6_a35; 
Det_3x3_LU_matlab(2) = (det([data_in(3:end,2), data_in(3:end,4:5)])); % 2 4 5 

DetM_3x3_n_13 = Mult_DetM_2x2_n_9_a32 - Mult_DetM_2x2_n_7_a33 + Mult_DetM_2x2_n_5_a35;
Det_3x3_LU_matlab(3) = (det([data_in(3:end,2), data_in(3:end,3), data_in(3:end,5)])); % 2 3 5

DetM_3x3_n_14 = Mult_DetM_2x2_n_8_a32 - Mult_DetM_2x2_n_6_a33 + Mult_DetM_2x2_n_5_a34;
Det_3x3_LU_matlab(4) = (det([data_in(3:end,2), data_in(3:end,3), data_in(3:end,4)])); 

% 2 3x3
DetM_3x3_n_22 = Mult_DetM_2x2_n_10_a31 - Mult_DetM_2x2_n_4_a34 + Mult_DetM_2x2_n_3_a35; 
Det_3x3_LU_matlab(5) = (det([data_in(3:end,1), data_in(3:end,4), data_in(3:end,5)])); 

DetM_3x3_n_23 = Mult_DetM_2x2_n_9_a31 - Mult_DetM_2x2_n_4_a33 + Mult_DetM_2x2_n_2_a35;
Det_3x3_LU_matlab(6) = (det([data_in(3:end,1), data_in(3:end,3), data_in(3:end,5)])); 

DetM_3x3_n_24 = Mult_DetM_2x2_n_8_a31 - Mult_DetM_2x2_n_3_a33 + Mult_DetM_2x2_n_2_a34;
Det_3x3_LU_matlab(8) = (det([data_in(3:end,1), data_in(3:end,2), data_in(3:end,5)])); 

% 3 3x3
DetM_3x3_n_33 = Mult_DetM_2x2_n_7_a31 - Mult_DetM_2x2_n_4_a32 + Mult_DetM_2x2_n_1_a35;
Det_3x3_LU_matlab(9) = (det([data_in(3:end,1), data_in(3:end,2), data_in(3:end,4)])); 

DetM_3x3_n_34 = Mult_DetM_2x2_n_6_a31 - Mult_DetM_2x2_n_3_a32 + Mult_DetM_2x2_n_1_a34;
Det_3x3_LU_matlab(7) = (det([data_in(3:end,1), data_in(3:end,3), data_in(3:end,4)])); 

% 4 3x3
DetM_3x3_n_44 = Mult_DetM_2x2_n_5_a31 - Mult_DetM_2x2_n_2_a32 + Mult_DetM_2x2_n_1_a33;
Det_3x3_LU_matlab(10) = (det([data_in(3:end,1), data_in(3:end,2), data_in(3:end,3)])); 

for n = 1:sim_options.num_det2x2
    if Det_3x3_LU_matlab(n) < 0
	    Det_3x3_LU_matlab(n) = Det_3x3_LU_matlab(n) * -1;
    end
end

%% Сумматоры определителя 3х3
[DetM_3x3_int_sum1(1,1), DetM_3x3_int_sum1_overflow(1,1), s.DetM_3x3_int_pre_sum_array(1,1), s.DetM_3x3_int_pre_sum_width_total(1,1)] = ...
    adder(Mult_DetM_3x3_int(1,1), -Mult_DetM_3x3_int(2,1), sim_options.type_3x3_det, sim_options.width_hilbert);

[DetM_3x3_int(1,1), DetM_3x3_int_overflow(1,1), s.DetM_3x3_int_sum_array(1,1), s.DetM_3x3_int_sum_width_total(1,1)] = ...
    adder(Mult_DetM_3x3_int(5,1), DetM_3x3_int_sum1(1,1), sim_options.type_3x3_det, sim_options.width_hilbert);
	
[DetM_3x3_int_sum1(2,1), DetM_3x3_int_sum1_overflow(2,1), s.DetM_3x3_int_pre_sum_array(2,1), s.DetM_3x3_int_pre_sum_width_total(2,1)] = ...
    adder(Mult_DetM_3x3_int(1,2), -Mult_DetM_3x3_int(3,1), sim_options.type_3x3_det, sim_options.width_hilbert);
[DetM_3x3_int(2,1), DetM_3x3_int_overflow(2,1), s.DetM_3x3_int_sum_array(2,1), s.DetM_3x3_int_sum_width_total(2,1)] = ...
    adder(Mult_DetM_3x3_int(6,1), DetM_3x3_int_sum1(2,1), sim_options.type_3x3_det, sim_options.width_hilbert);

[DetM_3x3_int_sum1(3,1), DetM_3x3_int_sum1_overflow(3,1), s.DetM_3x3_int_pre_sum_array(3,1), s.DetM_3x3_int_pre_sum_width_total(3,1)] = ...
    adder(Mult_DetM_3x3_int(1,3), -Mult_DetM_3x3_int(4,1), sim_options.type_3x3_det, sim_options.width_hilbert);
[DetM_3x3_int(3,1), DetM_3x3_int_overflow(3,1), s.DetM_3x3_int_sum_array(3,1), s.DetM_3x3_int_sum_width_total(3,1)] = ...
    adder(Mult_DetM_3x3_int(7,1), DetM_3x3_int_sum1(3,1), sim_options.type_3x3_det, sim_options.width_hilbert);

[DetM_3x3_int_sum1(4,1), DetM_3x3_int_sum1_overflow(4,1), s.DetM_3x3_int_pre_sum_array(4,1), s.DetM_3x3_int_pre_sum_width_total(4,1)] = ...
    adder(Mult_DetM_3x3_int(2,2), -Mult_DetM_3x3_int(3,2), sim_options.type_3x3_det, sim_options.width_hilbert);
[DetM_3x3_int(4,1), DetM_3x3_int_overflow(4,1), s.DetM_3x3_int_sum_array(4,1), s.DetM_3x3_int_sum_width_total(4,1)] = ...
    adder(Mult_DetM_3x3_int(8,1), DetM_3x3_int_sum1(4,1), sim_options.type_3x3_det, sim_options.width_hilbert);

[DetM_3x3_int_sum1(5,1), DetM_3x3_int_sum1_overflow(5,1), s.DetM_3x3_int_pre_sum_array(5,1), s.DetM_3x3_int_pre_sum_width_total(5,1)] = ...
    adder(Mult_DetM_3x3_int(2,3), -Mult_DetM_3x3_int(4,2), sim_options.type_3x3_det, sim_options.width_hilbert);
[DetM_3x3_int(5,1), DetM_3x3_int_overflow(5,1), s.DetM_3x3_int_sum_array(5,1), s.DetM_3x3_int_sum_width_total(5,1)] = ...
    adder(Mult_DetM_3x3_int(9,1), DetM_3x3_int_sum1(5,1), sim_options.type_3x3_det, sim_options.width_hilbert);

[DetM_3x3_int_sum1(6,1), DetM_3x3_int_sum1_overflow(6,1), s.DetM_3x3_int_pre_sum_array(6,1), s.DetM_3x3_int_pre_sum_width_total(6,1)] = ...
    adder(Mult_DetM_3x3_int(3,3), -Mult_DetM_3x3_int(4,3), sim_options.type_3x3_det, sim_options.width_hilbert);
[DetM_3x3_int(6,1), DetM_3x3_int_overflow(6,1), s.DetM_3x3_int_sum_array(6,1), s.DetM_3x3_int_sum_width_total(6,1)] = ...
    adder(Mult_DetM_3x3_int(10,1), DetM_3x3_int_sum1(6,1), sim_options.type_3x3_det, sim_options.width_hilbert);

[DetM_3x3_int_sum1(7,1), DetM_3x3_int_sum1_overflow(7,1), s.DetM_3x3_int_pre_sum_array(7,1), s.DetM_3x3_int_pre_sum_width_total(7,1)] = ...
    adder(Mult_DetM_3x3_int(5,2), -Mult_DetM_3x3_int(6,2), sim_options.type_3x3_det, sim_options.width_hilbert);
[DetM_3x3_int(7,1), DetM_3x3_int_overflow(7,1), s.DetM_3x3_int_sum_array(7,1), s.DetM_3x3_int_sum_width_total(7,1)] = ...
    adder(Mult_DetM_3x3_int(8,2), DetM_3x3_int_sum1(7,1), sim_options.type_3x3_det, sim_options.width_hilbert);

[DetM_3x3_int_sum1(8,1), DetM_3x3_int_sum1_overflow(8,1), s.DetM_3x3_int_pre_sum_array(8,1), s.DetM_3x3_int_pre_sum_width_total(8,1)] = ...
    adder(Mult_DetM_3x3_int(5,3), -Mult_DetM_3x3_int(7,2), sim_options.type_3x3_det, sim_options.width_hilbert);
[DetM_3x3_int(8,1), DetM_3x3_int_overflow(8,1), s.DetM_3x3_int_sum_array(8,1), s.DetM_3x3_int_sum_width_total(8,1)] = ...
    adder(Mult_DetM_3x3_int(9,2), DetM_3x3_int_sum1(8,1), sim_options.type_3x3_det, sim_options.width_hilbert);

[DetM_3x3_int_sum1(9,1), DetM_3x3_int_sum1_overflow(9,1), s.DetM_3x3_int_pre_sum_array(9,1), s.DetM_3x3_int_pre_sum_width_total(9,1)] = ...
    adder(Mult_DetM_3x3_int(6,3), -Mult_DetM_3x3_int(7,3), sim_options.type_3x3_det, sim_options.width_hilbert);
[DetM_3x3_int(9,1), DetM_3x3_int_overflow(9,1), s.DetM_3x3_int_sum_array(9,1), s.DetM_3x3_int_sum_width_total(9,1)] = ...
    adder(Mult_DetM_3x3_int(10,2), DetM_3x3_int_sum1(9,1), sim_options.type_3x3_det, sim_options.width_hilbert);

[DetM_3x3_int_sum1(10,1), DetM_3x3_int_sum1_overflow(10,1), s.DetM_3x3_int_pre_sum_array(10,1), s.DetM_3x3_int_pre_sum_width_total(10,1)] = ...
    adder(Mult_DetM_3x3_int(8,3), -Mult_DetM_3x3_int(9,3), sim_options.type_3x3_det, sim_options.width_hilbert);
[DetM_3x3_int(10,1), DetM_3x3_int_overflow(10,1), s.DetM_3x3_int_sum_array(10,1), s.DetM_3x3_int_sum_width_total(10,1)] = ...
    adder(Mult_DetM_3x3_int(10,3), DetM_3x3_int_sum1(10,1), sim_options.type_3x3_det, sim_options.width_hilbert);


for i = 1:sim_options.num_det2x2
    % Проверка переполнения пресумматора
    if (DetM_3x3_int_sum1_overflow(i,1) == 1)
		disp(['Sum det3x3 overflow in ', num2str(i)]);
        disp({DetM_3x3_int_sum1(i,1)});
    end
	% Проверка выходной разрядности пресумматора
    if (s.DetM_3x3_int_pre_sum_width_total(1,1) > sim_options.width_hilbert-1)
		disp(['Sum total width in function det3x3 higher than 64 in ', num2str(i)]);
        disp({DetM_3x3_int_sum1(i,1)});
    end
    % Наложение маски на пресумматор
    if sim_options.enable_mask == true
        c1 = bitmask(DetM_3x3_int_sum1(i,1), sim_options.type_3x3_det, det_max_width(i,4));
        if c1 ~= DetM_3x3_int_sum1(i,1)
            disp(['Bit mask error pre_sum det 3x3 in ', num2str(i)]);
            disp(DetM_3x3_int_sum1_width_total(i));
            disp({c1, DetM_3x3_int_sum1(i,1)});
            disp({sim_options.freq, sim_options.SNR});
        end
        DetM_3x3_int_sum1(i,1) = c1;
    end

    %% Проверка переполнения сумматора
    if (DetM_3x3_int_overflow(i,1) == 1)
		disp(['Sum det3x3 overflow in ', num2str(i)]);
        disp([DetM_3x3_int(i,1)]);
    end
	% Проверка выходной разрядности сумматора
    if (s.DetM_3x3_int_sum_width_total(i,1) > sim_options.width_hilbert-1)
		disp(['Sum total width in function det3x3 higher than 64 in ', num2str(i)]);
        disp([DetM_3x3_int(i,1)]);
    end
    % Наложение маски на сумматор
    if sim_options.enable_mask == true
        c1 = bitmask(DetM_3x3_int(i,1), sim_options.type_3x3_det, det_max_width(i,5));
        if c1 ~= DetM_3x3_int(i,1)
            disp(['Bit mask error pre_sum det 3x3 in ', num2str(i)]);
            % c2 = bitmask(sum(i+1,n), sim_options.int_size, width_sum(i+1));
            disp(s.DetM_3x3_int_sum_width_total(i));
            disp([c1, DetM_3x3_int(i,1)]);
            disp([sim_options.freq, sim_options.SNR]);
        end
        DetM_3x3_int(i,1) = c1;
    end
end

                                                                            %% 4x4
% double                                                            
var_mult_a_det11 = data_in(2,2) * DetM_3x3_n_11; 
var_mult_a_det12 = data_in(2,3) * DetM_3x3_n_12;
var_mult_a_det13 = data_in(2,4) * DetM_3x3_n_13;
var_mult_a_det14 = data_in(2,5) * DetM_3x3_n_14;

DetM_4x4_n_1 = var_mult_a_det11 - var_mult_a_det12  + var_mult_a_det13 - var_mult_a_det14;
Det_4x4_LU_matlab(1) = det([data_in(2:end,2), data_in(2:end,3), data_in(2:end,4), data_in(2:end,5)]); 

var_mult_a_det21 = data_in(2,1) * DetM_3x3_n_11;
var_mult_a_det22 = data_in(2,3) * DetM_3x3_n_22;
var_mult_a_det23 = data_in(2,4) * DetM_3x3_n_23;
var_mult_a_det24 = data_in(2,5) * DetM_3x3_n_24;

DetM_4x4_n_2 = var_mult_a_det21 - var_mult_a_det22 + var_mult_a_det23 - var_mult_a_det24;
Det_4x4_LU_matlab(2) = det([data_in(2:end,1), data_in(2:end,3), data_in(2:end,4), data_in(2:end,5)]); 

var_mult_a_det31 = data_in(2,1) * DetM_3x3_n_12;
var_mult_a_det32 = data_in(2,2) * DetM_3x3_n_22;
var_mult_a_det33 = data_in(2,4) * DetM_3x3_n_33;
var_mult_a_det34 = data_in(2,5) * DetM_3x3_n_34;

DetM_4x4_n_3 = var_mult_a_det31 - var_mult_a_det32 + var_mult_a_det33 - var_mult_a_det34;
Det_4x4_LU_matlab(3) = det([data_in(2:end,1), data_in(2:end,2), data_in(2:end,4), data_in(2:end,5)]); 

var_mult_a_det41 = data_in(2,1) * DetM_3x3_n_13;
var_mult_a_det42 = data_in(2,2) * DetM_3x3_n_23;
var_mult_a_det43 = data_in(2,3) * DetM_3x3_n_33;
var_mult_a_det44 = data_in(2,5) * DetM_3x3_n_44;

DetM_4x4_n_4 = var_mult_a_det41 - var_mult_a_det42 + var_mult_a_det43 - var_mult_a_det44;
Det_4x4_LU_matlab(4) = det([data_in(2:end,1), data_in(2:end,2), data_in(2:end,3), data_in(2:end,5)]); 

var_mult_a_det51 = data_in(2,1) * DetM_3x3_n_14;
var_mult_a_det52 = data_in(2,2) * DetM_3x3_n_24;
var_mult_a_det53 = data_in(2,3) * DetM_3x3_n_34;
var_mult_a_det54 = data_in(2,4) * DetM_3x3_n_44;

DetM_4x4_n_5 = var_mult_a_det51 - var_mult_a_det52 + var_mult_a_det53 - var_mult_a_det54;
Det_4x4_LU_matlab(5) = det([data_in(2:end,1), data_in(2:end,2), data_in(2:end,3), data_in(2:end,4)]); 

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
        [var_mult_det4x4_int(j,i), var_mult_det4x4_int_overflow(j,i), var_mult_det4x4_abs(j,i), var_mult_det4x4_width_total(j,i)] = mult(cast(data_in_int(2,index1(j,i)), sim_options.type_4x4_det), ...
            cast(DetM_3x3_int(index_DetM_3x3(j,i)), sim_options.type_4x4_det), sim_options.type_4x4_det, sim_options.width_hilbert);

        % Проверка переполнения умножителя
        if (var_mult_det4x4_int_overflow(j,i) == 1)
		    disp(['Mult det4x4 overflow in ', num2str(j), num2str(i)]);
            disp({var_mult_det4x4_int(j,i)});
        end
	    % Проверка выходной разрядности умножителя
        if (var_mult_det4x4_width_total(j,i) > sim_options.width_hilbert-1)
		    disp(['Mult total width in function det4x4 higher than 64 in ', num2str(j), num2str(i)]);
            disp({var_mult_det4x4_int(j,i)});
        end       
    end
    s.DetM_4x4_int_mult_array(kk+1:kk+4) = var_mult_det4x4_abs(j,:);
    s.DetM_4x4_int_mult_width_total(kk+1:kk+4) = var_mult_det4x4_width_total(j,:);
    DetM_4x4_int_mult_array(kk+1:kk+4) = var_mult_det4x4_abs(j,:);
    kk = kk + 4;
end

% Наложение маски на умножитель
for j = 1:sim_options.num_det2x2*2
    if sim_options.enable_mask == true
        c1 = bitmask(DetM_4x4_int_mult_array(j), sim_options.type_3x3_det, det_max_width(j,6));
        if c1 ~= DetM_4x4_int_mult_array(j)
            disp(['Bit mask error pre_sum det 3x3 in ', num2str(j)]);
            % c2 = bitmask(sum(i+1,n), sim_options.int_size, width_sum(i+1));
            disp([c1, DetM_4x4_int_mult_array(j)]);
            disp([sim_options.freq, sim_options.SNR]);
        end
        DetM_4x4_int_mult_array(j) = c1;
    end 
end

%% Сумматоры определителя 4х4
kk = 1;
for i = 1:5
    for j = 1:2
        if (mod(j,2) == 1)
            [DetM_4x4_int_sum(i,j), DetM_4x4_int_sum1_overflow(i,j), DetM_4x4_int_sum1_abs(i,j), DetM_4x4_int_sum1_width_total(i,j)] = adder(var_mult_det4x4_int(i,j), ...
                -var_mult_det4x4_int(i,j+1),  sim_options.type_4x4_det, sim_options.width_hilbert);
        
            % Проверка переполнения сумматора
            if (DetM_4x4_int_sum1_overflow(i,j) == 1)
		        disp(['Sum1 det4x4 overflow in ', num2str(i), num2str(j)]);
                disp({DetM_4x4_int_sum(i,j)});
            end
	        % Проверка выходной разрядности сумматора
            if (DetM_4x4_int_sum1_width_total(i,j) > sim_options.width_hilbert-1)
		        disp(['Sum1 total width in function det4x4 higher than 64 in ', num2str(i), num2str(j)]);
                disp({DetM_4x4_int_sum(i,1)});
            end
        else
            [DetM_4x4_int_sum(i,j), DetM_4x4_int_sum1_overflow(i,j), DetM_4x4_int_sum1_abs(i,j), DetM_4x4_int_sum1_width_total(i,j)] = adder(var_mult_det4x4_int(i,j+1), ...
                -var_mult_det4x4_int(i,j+2),  sim_options.type_4x4_det, sim_options.width_hilbert);
         
            % Проверка переполнения сумматора
            if (DetM_4x4_int_sum1_overflow(i,j) == 1)
		        disp(['Sum2 det4x4 overflow in ', num2str(i), num2str(j)]);
                disp({DetM_4x4_int_sum(i,j)});
            end
	        % Проверка выходной разрядности сумматора
            if (DetM_4x4_int_sum1_width_total(i,j) > sim_options.width_hilbert-1)
		        disp(['Sum2 total width in function det4x4 higher than 64 in ', num2str(i), num2str(j)]);
                disp({DetM_4x4_int_sum(i,j)});
            end
        end
        s.DetM_4x4_int_pre_sum_array(kk) = DetM_4x4_int_sum1_abs(i,j);
        s.DetM_4x4_int_pre_sum_width_total(kk) = DetM_4x4_int_sum1_width_total(i,j);
        DetM_4x4_int_pre_sum_array(kk) = DetM_4x4_int_sum(i,j);
        
        % Наложение маски на сумматор
        if sim_options.enable_mask == true
            c1 = bitmask(DetM_4x4_int_pre_sum_array(kk), sim_options.type_4x4_det, det_max_width(kk,7));
            if c1 ~= DetM_4x4_int_pre_sum_array(kk)
                disp(['Bit mask error pre_sum det 4x4 in ', num2str(kk)]);
                disp({c1, DetM_4x4_int_pre_sum_array(kk)});
                disp({sim_options.freq, sim_options.SNR});
            end
            DetM_4x4_int_pre_sum_array(i,j) = c1;
        end
        kk = kk + 1;
    end


    % Сумматор определителя 4х4
    [DetM_4x4_int(i), DetM_4x4_int_overflow(i), s.DetM_4x4_int_sum_array(i), s.DetM_4x4_int_sum_array_width_total(i)] = adder(DetM_4x4_int_sum(i,1), ...
        DetM_4x4_int_sum(i,2), sim_options.type_4x4_det, sim_options.width_hilbert);

    % Проверка переполнения сумматора
    if (DetM_4x4_int_overflow(i) == 1)
	    disp(['Sum det4x4 overflow in ', num2str(i)]);
        disp([DetM_4x4_int(i)]);
    end
	% Проверка выходной разрядности сумматора
    if (s.DetM_4x4_int_sum_array_width_total(i) > sim_options.width_hilbert-1)
	    disp(['Sum total width in function det4x4 higher than 64 in ', num2str(i)]);
        disp([DetM_4x4_int(i)]);
    end

    % Наложение маски на сумматор
    if sim_options.enable_mask == true
        c1 = bitmask(DetM_4x4_int(i), sim_options.type_4x4_det, det_max_width(i,8));
        if c1 ~= DetM_4x4_int(i)
            disp(['Bit mask error sum det 4x4 in ', num2str(i)]);
            disp([c1, DetM_4x4_int(i)]);
            disp([sim_options.freq, sim_options.SNR]);
        end
        DetM_4x4_int(i) = c1;
    end
end

                                                                            %% 5x5
	% double                                                                    
	var1 = data_in(1,1) * DetM_4x4_n_1;
	var2 = data_in(1,2) * DetM_4x4_n_2;
	var3 = data_in(1,3) * DetM_4x4_n_3;
	var4 = data_in(1,4) * DetM_4x4_n_4;
	var5 = data_in(1,5) * DetM_4x4_n_5;
	DetM_5x5 = var1 - var2 + var3 - var4 + var5;

	%% Умножители матрицы 5х5
	for i = 1:sim_options.Size_matrix

		[var_int_det5x5(i,:), var_int_det5x5_overflow(i,:), s.DetM_5x5_int_mult_array(i,:), s.DetM_5x5_int_mult_array_width_total(i,:)] = mult(cast(data_in_int(1,i),sim_options.type_5x5_det), ...
			DetM_4x4_int(i), sim_options.type_5x5_det, sim_options.width_hilbert);
	
		% Проверка переполнения сумматора
		if (var_int_det5x5_overflow(i) == 1)
			disp(['Mult det5x5 overflow in ', num2str(i)]);
			disp([var_int_det5x5(i)]);
		end
		% Проверка выходной разрядности сумматора
		if (s.DetM_5x5_int_mult_array_width_total(i) > sim_options.width_hilbert-1)
			disp(['Mult total width in function det5x5 higher than 64 in ', num2str(i)]);
			disp([var_int_det5x5(i)]);
		end

		% Наложение маски на сумматор
		if sim_options.enable_mask == true
			c1 = bitmask(var_int_det5x5(i), sim_options.type_5x5_det, det_max_width(i,9));
			if c1 ~= var_int_det5x5(i)
				disp(['Bit mask error mult det 5x5 in ', num2str(i)]);
				% c2 = bitmask(sum(i+1,n), sim_options.int_size, width_sum(i+1));
				% disp(DetM_4x4_int_sum(i,j));
				disp([c1, var_int_det5x5(i)]);
				disp([sim_options.freq, sim_options.SNR]);
			end
			var_int_det5x5(i) = c1;
        end
	end

	%% Сумматоры1 определителя 5х5
	for i = 1:2
		if (mod(i,2) == 1)
			[DetM_5x5_int_sum1(i), DetM_5x5_int_sum1_overflow(i), s.DetM_5x5_int_pre_sum1_array(i), s.DetM_5x5_int_pre_sum1_array_width_total(i)] = adder(var_int_det5x5(i), ...
				-var_int_det5x5(i+1), sim_options.type_5x5_det, sim_options.width_hilbert);
			% Проверка переполнения сумматора
			if (DetM_5x5_int_sum1_overflow(i) == 1)
				disp(['Sum1 det5x5 overflow in ', num2str(i)]);
				disp([DetM_5x5_int_sum1_overflow(i)]);
			end
			% Проверка выходной разрядности сумматора
			if (s.DetM_5x5_int_pre_sum1_array_width_total(i) > sim_options.width_hilbert-1)
				disp(['Sum1 total width in function det5x5 higher than 64 in ', num2str(i)]);
				disp([DetM_5x5_int_sum1_width_total(i)]);
            end
		else
			[DetM_5x5_int_sum1(i), DetM_5x5_int_sum1_overflow(i), s.DetM_5x5_int_pre_sum1_array(i), s.DetM_5x5_int_pre_sum1_array_width_total(i)] = adder(var_int_det5x5(i+1), ...
				-var_int_det5x5(i+2), sim_options.type_5x5_det, sim_options.width_hilbert);
			% Проверка переполнения сумматора
			if (DetM_5x5_int_sum1_overflow(i) == 1)
				disp(['Sum1 det5x5 overflow in ', num2str(i)]);
				disp([DetM_5x5_int_sum1(i)]);
			end
			% Проверка выходной разрядности сумматора
			if (s.DetM_5x5_int_pre_sum1_array_width_total(i) > sim_options.width_hilbert-1)
				disp(['Sum1 total width in function det5x5 higher than 64 in ', num2str(i)]);
				disp([DetM_5x5_int_sum1(i)]);
            end
        end

        % Наложение маски на сумматор
		if sim_options.enable_mask == true
		    c1 = bitmask(DetM_5x5_int_sum1(i), sim_options.type_5x5_det, det_max_width(i,10));
			if c1 ~= DetM_5x5_int_sum1(i)
			    disp(['Bit mask error sum det 5x5 in ', num2str(i)]);
				disp([c1, DetM_5x5_int_sum1(i)]);
				disp([sim_options.freq, sim_options.SNR]);
            end
			DetM_5x5_int_sum1(i) = c1;
        end
	end

	%% Сумматор2 определителя 5х5
	[DetM_5x5_int_sum3, DetM_5x5_int_sum3_overflow, s.DetM_5x5_int_sum3_abs, s.DetM_5x5_int_sum3_width_total] = adder(DetM_5x5_int_sum1(1), ...
		DetM_5x5_int_sum1(2), sim_options.type_5x5_det, sim_options.width_hilbert);
	
	% Проверка переполнения сумматора
	if (DetM_5x5_int_sum3_overflow == 1)
		disp('Sum3 det5x5 overflow');
		disp({DetM_5x5_int_sum3});
	end
	% Проверка выходной разрядности сумматора
	if (s.DetM_5x5_int_sum3_width_total > sim_options.width_hilbert-1)
		disp('Sum3 total width in function det5x5 higher than 64');
		disp({DetM_5x5_int_sum3});
	end

	% Наложение маски на сумматор
	if sim_options.enable_mask == true
		c1 = bitmask(DetM_5x5_int_sum3, sim_options.type_5x5_det, det_max_width(1,11));
		if c1 ~= DetM_5x5_int_sum3
			disp('Bit mask error sum3 det 5x5');
			disp([c1, DetM_5x5_int_sum3]);
			disp([sim_options.freq, sim_options.SNR]);
		end
		DetM_5x5_int_sum3 = c1;
	end

	%% Сумматор3 определителя 5х5
	[DetM_5x5_int, DetM_5x5_int_overflow, s.DetM_5x5_int_abs, s.DetM_5x5_int_width_total] = adder(DetM_5x5_int_sum3,  var_int_det5x5(5), ...
		sim_options.type_5x5_det, sim_options.width_hilbert);

	% Проверка переполнения сумматора
	if (DetM_5x5_int_overflow == 1)
		disp('Sum3 det5x5 overflow');
		disp([DetM_5x5_int]);
	end
	% Проверка выходной разрядности сумматора
	if (s.DetM_5x5_int_width_total > sim_options.width_hilbert-1)
		disp('Sum3 total width in function det5x5 higher than 64');
		disp([DetM_5x5_int]);
	end
	% Наложение маски на сумматор
	if sim_options.enable_mask == true
		c1 = bitmask(DetM_5x5_int, sim_options.type_5x5_det, det_max_width(1,11));
		if c1 ~= DetM_5x5_int
			disp('Bit mask error sum det 5x5');
			disp({c1, DetM_5x5_int});
			disp({sim_options.freq, sim_options.SNR});
		end
		DetM_5x5_int = c1;
	end

    %% Если определитель равен нулю, то присваиваем значение 1. Обязательно!
    if (DetM_5x5_int == 0)
        DetM_5x5_int = 1;
        disp('Determinant int x3 equal 0');
    end
    
    if (DetM_5x5 == 0)
        DetM_5x5 = 1;
        disp('Determinant double x3 equal 0');
    end



	%% double
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