
function [DetM_5x5, DetM_5x5_int, DetM_2x2, DetM_2x2_int, Det2x2_mult1_abs, Det2x2_mult2_abs, Det2x2_sum_abs, ...
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
] = determinate(data_in, data_in_int, int_size, width)

a = data_in;
a_int = data_in_int;

int_size_double = "double";
width_double = 80;
%% Находим матрицы 2x2
s = struct;
e = 0;
t=4;
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

DetM_2x2 = zeros(10,1);
DetM_2x2_int = cast(zeros(10,1), int_size);

Mult_DetM_2x2_n_1_int = cast(zeros(3,1),  int_size_double);
Mult_DetM_2x2_n_2_int = cast(zeros(3,1),  int_size_double);
Mult_DetM_2x2_n_3_int = cast(zeros(3,1),  int_size_double);
Mult_DetM_2x2_n_4_int = cast(zeros(3,1),  int_size_double);
Mult_DetM_2x2_n_5_int = cast(zeros(3,1),  int_size_double);
Mult_DetM_2x2_n_6_int = cast(zeros(3,1),  int_size_double);
Mult_DetM_2x2_n_7_int = cast(zeros(3,1),  int_size_double);
Mult_DetM_2x2_n_8_int = cast(zeros(3,1),  int_size_double);
Mult_DetM_2x2_n_9_int = cast(zeros(3,1),  int_size_double);
Mult_DetM_2x2_n_10_int = cast(zeros(3,1), int_size_double);

Mult_DetM_2x2_n_1_overflow = cast(zeros(3,1),  int_size);
Mult_DetM_2x2_n_2_overflow = cast(zeros(3,1),  int_size);
Mult_DetM_2x2_n_3_overflow = cast(zeros(3,1),  int_size);
Mult_DetM_2x2_n_4_overflow = cast(zeros(3,1),  int_size);
Mult_DetM_2x2_n_5_overflow = cast(zeros(3,1),  int_size);
Mult_DetM_2x2_n_6_overflow = cast(zeros(3,1),  int_size);
Mult_DetM_2x2_n_7_overflow = cast(zeros(3,1),  int_size);
Mult_DetM_2x2_n_8_overflow = cast(zeros(3,1),  int_size);
Mult_DetM_2x2_n_9_overflow = cast(zeros(3,1),  int_size);
Mult_DetM_2x2_n_10_overflow = cast(zeros(3,1), int_size);

Mult_DetM_2x2_n_1_abs = cast(zeros(3,1),  int_size_double);
Mult_DetM_2x2_n_2_abs = cast(zeros(3,1),  int_size_double);
Mult_DetM_2x2_n_3_abs = cast(zeros(3,1),  int_size_double);
Mult_DetM_2x2_n_4_abs = cast(zeros(3,1),  int_size_double);
Mult_DetM_2x2_n_5_abs = cast(zeros(3,1),  int_size_double);
Mult_DetM_2x2_n_6_abs = cast(zeros(3,1),  int_size_double);
Mult_DetM_2x2_n_7_abs = cast(zeros(3,1),  int_size_double);
Mult_DetM_2x2_n_8_abs = cast(zeros(3,1),  int_size_double);
Mult_DetM_2x2_n_9_abs = cast(zeros(3,1),  int_size_double);
Mult_DetM_2x2_n_10_abs = cast(zeros(3,1), int_size_double);

Mult_DetM_2x2_n_1_total_width = cast(zeros(3,1),  int_size);
Mult_DetM_2x2_n_2_total_width = cast(zeros(3,1),  int_size);
Mult_DetM_2x2_n_3_total_width = cast(zeros(3,1),  int_size);
Mult_DetM_2x2_n_4_total_width = cast(zeros(3,1),  int_size);
Mult_DetM_2x2_n_5_total_width = cast(zeros(3,1),  int_size);
Mult_DetM_2x2_n_6_total_width = cast(zeros(3,1),  int_size);
Mult_DetM_2x2_n_7_total_width = cast(zeros(3,1),  int_size);
Mult_DetM_2x2_n_8_total_width = cast(zeros(3,1),  int_size);
Mult_DetM_2x2_n_9_total_width = cast(zeros(3,1),  int_size);
Mult_DetM_2x2_n_10_total_width = cast(zeros(3,1), int_size);

Det2x2_mult1_abs = cast(zeros(10,1), int_size);
Det2x2_mult2_abs = cast(zeros(10,1), int_size);
Det2x2_sum_abs = cast(zeros(10,1), 	 int_size);


%% Находим определители матриц 2x2
for i = 1:10
    [DetM_2x2(i), DetM_2x2_int(i), Det2x2_mult1_abs(i), Det2x2_mult2_abs(i), Det2x2_sum_abs(i), mult1_overflow(i), mult2_overflow(i), width_total_mult1(i), ...
     width_total_mult2(i), sum_overflow(i), width_total_sum(i)]  = det_2x2(s.a{i}, s.a_int{i}, int_size, width); 

    %% Проверка переполнения умножителя
    if (mult1_overflow(i) == 1 || mult2_overflow(i) == 1)
		disp('Mult in function det2x2 overflow');
        disp(width);
        disp({s.a_int{i}});
    end
    % 
	%% Проверка выходной разрядности умножителя
    if (width_total_mult1(i) > width || width_total_mult2(i) > width)
		disp('Mult total width in function det2x2 higher than 64');
        disp(width);
        disp({s.a_int{i}});
    end
    % 
	%% Проверка переполнения сумматора
    if (sum_overflow(i) == 1)
		disp('Adder in function det2x2 overflow');
        disp(width);
        disp({s.a_int{i}});
    end
    % 
	%% Проверка выходной разрядности сумматора
    if (width_total_sum(i) > width)
		disp('Sum total width in function det2x2 higher than 64');
        disp(width);
        disp({width_total_sum(i)});
    end
	
end

%% double

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
    [Mult_DetM_2x2_n_7_int(i), Mult_DetM_2x2_n_7_overflow(i), Mult_DetM_2x2_n_7_abs(i), Mult_DetM_2x2_n_7_total_width(i)] = mult(double(a_int(3,index1(i))), double(DetM_2x2_int(7)), int_size_double, width_double);
    %% Проверка переполнения умножителя
    if (Mult_DetM_2x2_n_7_overflow(i) == 1)
		disp('Mult in function det3x3 7 overflow');
        disp({Mult_DetM_2x2_n_7_int(i)});
    end
	%% Проверка выходной разрядности умножителя
    if (Mult_DetM_2x2_n_7_total_width(i) > width_double-1)
		disp('Mult total width in function det3x3 7 higher than 64');
        disp({Mult_DetM_2x2_n_7_int(i)});
    end
    %%
    [Mult_DetM_2x2_n_9_int(i), Mult_DetM_2x2_n_9_overflow(i), Mult_DetM_2x2_n_9_abs(i), Mult_DetM_2x2_n_9_total_width(i)] = mult(double(a_int(3,index2(i))), double(DetM_2x2_int(9)), int_size_double, width_double);
    %% Проверка переполнения умножителя
    if (Mult_DetM_2x2_n_9_overflow(i) == 1)
		disp('Mult in function det3x3 9 overflow');
        disp({Mult_DetM_2x2_n_9_int(i)});
    end
	%% Проверка выходной разрядности умножителя
    if (Mult_DetM_2x2_n_9_total_width(i) > width_double-1)
		disp('Mult total width in function det3x3 9 higher than 64');
        disp({Mult_DetM_2x2_n_9_int(i)});
    end
	%%
    [Mult_DetM_2x2_n_8_int(i), Mult_DetM_2x2_n_8_overflow(i), Mult_DetM_2x2_n_8_abs(i), Mult_DetM_2x2_n_8_total_width(i)] = mult(double(a_int(3,index3(i))), double(DetM_2x2_int(8)), int_size_double, width_double);
    %% Проверка переполнения умножителя
    if (Mult_DetM_2x2_n_8_overflow(i) == 1)
		disp('Mult in function det3x3 8 overflow');
        disp(width);
        disp({Mult_DetM_2x2_n_8_int(i)});
    end
	%% Проверка выходной разрядности умножителя
    if (Mult_DetM_2x2_n_8_total_width(i) > width_double-1)
		disp('Mult total width in function det3x3 8 higher than 64');
        disp(width);
        disp({Mult_DetM_2x2_n_8_int(i)});
    end
    %%
    [Mult_DetM_2x2_n_6_int(i), Mult_DetM_2x2_n_6_overflow(i), Mult_DetM_2x2_n_6_abs(i), Mult_DetM_2x2_n_6_total_width(i)] = mult(double(a_int(3,index4(i))), double(DetM_2x2_int(6)), int_size_double, width_double);
    %% Проверка переполнения умножителя
    if (Mult_DetM_2x2_n_6_overflow(i) == 1)
		disp('Mult in function det3x3 6 overflow');
        disp(width);
        disp({Mult_DetM_2x2_n_6_int(i)});
    end
	%% Проверка выходной разрядности умножителя
    if (Mult_DetM_2x2_n_6_total_width(i) > width_double-1)
		disp('Mult total width in function det3x3 6 higher than 64');
        disp(width);
        disp({Mult_DetM_2x2_n_6_int(i)});
    end
	%%
    [Mult_DetM_2x2_n_5_int(i), Mult_DetM_2x2_n_5_overflow(i), Mult_DetM_2x2_n_5_abs(i), Mult_DetM_2x2_n_5_total_width(i)] = mult(double(a_int(3,index5(i))), double(DetM_2x2_int(5)), int_size_double, width_double);
    %% Проверка переполнения умножителя
    if (Mult_DetM_2x2_n_5_overflow(i) == 1)
		disp('Mult in function det3x3 5 overflow');
        disp(width);
        disp({Mult_DetM_2x2_n_5_int(i)});
    end
	%% Проверка выходной разрядности умножителя
    if (Mult_DetM_2x2_n_5_total_width(i) > width_double-1)
		disp('Mult total width in function det3x3 5 higher than 64');
        disp(width);
        disp({Mult_DetM_2x2_n_5_int(i)});
    end
    %%
    [Mult_DetM_2x2_n_4_int(i), Mult_DetM_2x2_n_4_overflow(i), Mult_DetM_2x2_n_4_abs(i), Mult_DetM_2x2_n_4_total_width(i)] = mult(double(a_int(3,index6(i))), double(DetM_2x2_int(4)), int_size_double, width_double);
    %% Проверка переполнения умножителя
    if (Mult_DetM_2x2_n_4_overflow(i) == 1)
		disp('Mult in function det3x3 4 overflow');
        disp(width);
        disp({Mult_DetM_2x2_n_4_int(i)});
    end
	%% Проверка выходной разрядности умножителя
    if (Mult_DetM_2x2_n_4_total_width(i) > width_double-1)
		disp('Mult total width in function det3x3 4 higher than 64');
        disp(width);
        disp({Mult_DetM_2x2_n_4_int(i)});
    end
    %%  
    [Mult_DetM_2x2_n_3_int(i), Mult_DetM_2x2_n_3_overflow(i), Mult_DetM_2x2_n_3_abs(i), Mult_DetM_2x2_n_3_total_width(i)] = mult(double(a_int(3,index7(i))), double(DetM_2x2_int(3)), int_size_double, width_double);
    %% Проверка переполнения умножителя
    if (Mult_DetM_2x2_n_3_overflow(i) == 1)
		disp('Mult in function det3x3 3 overflow');
        disp(width);
        disp({Mult_DetM_2x2_n_3_int(i)});
    end
	%% Проверка выходной разрядности умножителя
    if (Mult_DetM_2x2_n_3_total_width(i) > width_double-1)
		disp('Mult total width in function det3x3 3 higher than 64');
        disp(width);
        disp({Mult_DetM_2x2_n_3_int(i)});
    end
    %%
    [Mult_DetM_2x2_n_2_int(i), Mult_DetM_2x2_n_2_overflow(i), Mult_DetM_2x2_n_2_abs(i), Mult_DetM_2x2_n_2_total_width(i)] = mult(double(a_int(3,index8(i))), double(DetM_2x2_int(2)), int_size_double, width_double);
    %% Проверка переполнения умножителя
    if (Mult_DetM_2x2_n_2_overflow(i) == 1)
		disp('Mult in function det3x3 2 overflow');
        disp(width);
        disp({Mult_DetM_2x2_n_2_int(i)});
    end
	%% Проверка выходной разрядности умножителя
    if (Mult_DetM_2x2_n_2_total_width(i) > width_double-1)
		disp('Mult total width in function det3x3 2 higher than 64');
        disp(width);
        disp({Mult_DetM_2x2_n_2_int(i)});
    end
    %%
    [Mult_DetM_2x2_n_1_int(i), Mult_DetM_2x2_n_1_overflow(i), Mult_DetM_2x2_n_1_abs(i), Mult_DetM_2x2_n_1_total_width(i)] = mult(double(a_int(3,index9(i))), double(DetM_2x2_int(1)), int_size_double, width_double);
    %% Проверка переполнения умножителя
    if (Mult_DetM_2x2_n_1_overflow(i) == 1)
		disp('Mult in function det3x3 1 overflow');
        disp(width);
        disp({Mult_DetM_2x2_n_1_int(i)});
    end
	%% Проверка выходной разрядности умножителя
    if (Mult_DetM_2x2_n_1_total_width(i) > width_double-1)
		disp('Mult total width in function det3x3 1 higher than 64');
        disp(width);
        disp({Mult_DetM_2x2_n_1_int(i)});
    end
    %%
    [Mult_DetM_2x2_n_10_int(i), Mult_DetM_2x2_n_10_overflow(i), Mult_DetM_2x2_n_10_abs(i), Mult_DetM_2x2_n_10_total_width(i)] = mult(double(a_int(3,i)), double(DetM_2x2_int(10)), int_size_double, width_double);
    %% Проверка переполнения умножителя
    if (Mult_DetM_2x2_n_10_overflow(i) == 1)
		disp('Mult in function det3x3 10 overflow');
        disp(width);
        disp({Mult_DetM_2x2_n_10_int(i)});
    end
	%% Проверка выходной разрядности умножителя
    if (Mult_DetM_2x2_n_7_total_width(i) > width_double-1)
		disp('Mult total width in function det3x3 10 higher than 64');
        disp(width);
        disp({Mult_DetM_2x2_n_10_int(i)});
    end
end

%%
% 1 4x4
DetM_3x3_n_11 = Mult_DetM_2x2_n_10_a33 - Mult_DetM_2x2_n_9_a34 + Mult_DetM_2x2_n_8_a35; 
DetM_3x3_n_12 = Mult_DetM_2x2_n_10_a32 - Mult_DetM_2x2_n_7_a34 + Mult_DetM_2x2_n_6_a35; 
DetM_3x3_n_13 = Mult_DetM_2x2_n_9_a32 - Mult_DetM_2x2_n_7_a33 + Mult_DetM_2x2_n_5_a35;
DetM_3x3_n_14 = Mult_DetM_2x2_n_8_a32 - Mult_DetM_2x2_n_6_a33 + Mult_DetM_2x2_n_5_a34;
% 2 4x4
DetM_3x3_n_22 = Mult_DetM_2x2_n_10_a31 - Mult_DetM_2x2_n_4_a34 + Mult_DetM_2x2_n_3_a35; 
DetM_3x3_n_23 = Mult_DetM_2x2_n_9_a31 - Mult_DetM_2x2_n_4_a33 + Mult_DetM_2x2_n_2_a35;
DetM_3x3_n_24 = Mult_DetM_2x2_n_8_a31 - Mult_DetM_2x2_n_3_a33 + Mult_DetM_2x2_n_2_a34;
% 3 4x4
DetM_3x3_n_33 = Mult_DetM_2x2_n_7_a31 - Mult_DetM_2x2_n_4_a32 + Mult_DetM_2x2_n_1_a35;
DetM_3x3_n_34 = Mult_DetM_2x2_n_6_a31 - Mult_DetM_2x2_n_3_a32 + Mult_DetM_2x2_n_1_a34;
% 4 4x4
DetM_3x3_n_44 = Mult_DetM_2x2_n_5_a31 - Mult_DetM_2x2_n_2_a32 + Mult_DetM_2x2_n_1_a33;

%%
% 4x4 
[DetM_3x3_n_11_int_sum1, DetM_3x3_n_11_int_sum1_overflow, DetM_3x3_n_11_int_sum1_abs, DetM_3x3_n_11_int_sum1_width_total] = ...
    adder(Mult_DetM_2x2_n_10_int(3),  -Mult_DetM_2x2_n_9_int(3), int_size_double, width_double);
[DetM_3x3_n_11_int, DetM_3x3_n_11_int_overflow, DetM_3x3_n_11_int_abs, DetM_3x3_n_11_int_width_total] = ...
    adder(DetM_3x3_n_11_int_sum1,  Mult_DetM_2x2_n_8_int(3), int_size_double, width_double);
    %% Проверка переполнения сумматора
    if (DetM_3x3_n_11_int_overflow == 1)
		disp('Sum det3x3 11 overflow');
        disp(width);
        disp({DetM_3x3_n_11_int});
    end
	%% Проверка выходной разрядности сумматора
    if (DetM_3x3_n_11_int_width_total > width_double-1)
		disp('Sum total width in function det3x3 11 higher than 64');
        disp(width);
        disp({DetM_3x3_n_11_int});
    end
    %%
[DetM_3x3_n_12_int_sum1, DetM_3x3_n_12_int_sum1_overflow, DetM_3x3_n_12_int_sum1_abs, DetM_3x3_n_12_int_sum1_width_total] = ...
    adder(Mult_DetM_2x2_n_10_int(2),  -Mult_DetM_2x2_n_7_int(3), int_size_double, width_double);
[DetM_3x3_n_12_int, DetM_3x3_n_12_int_overflow, DetM_3x3_n_12_int_abs, DetM_3x3_n_12_int_width_total] = ...
    adder(DetM_3x3_n_12_int_sum1,  Mult_DetM_2x2_n_6_int(3), int_size_double, width_double);
    %% Проверка переполнения сумматора
    if (DetM_3x3_n_12_int_overflow == 1)
		disp('Sum det3x3 12 overflow');
        disp(width);
        disp({DetM_3x3_n_12_int});
    end
	%% Проверка выходной разрядности сумматора
    if (DetM_3x3_n_12_int_width_total > width_double-1)
		disp('Sum total width in function det3x3 12 higher than 64');
        disp(width);
        disp({DetM_3x3_n_12_int});
    end
    %%
[DetM_3x3_n_13_int_sum1, DetM_3x3_n_13_int_sum1_overflow, DetM_3x3_n_13_int_sum1_abs, DetM_3x3_n_13_int_sum1_width_total] = ...
    adder(Mult_DetM_2x2_n_9_int(2),  -Mult_DetM_2x2_n_7_int(2), int_size_double, width_double);
[DetM_3x3_n_13_int, DetM_3x3_n_13_int_overflow, DetM_3x3_n_13_int_abs, DetM_3x3_n_13_int_width_total] = ...
    adder(DetM_3x3_n_13_int_sum1,  Mult_DetM_2x2_n_5_int(3), int_size_double, width_double);
    %% Проверка переполнения сумматора
    if (DetM_3x3_n_13_int_overflow == 1)
		disp('Sum det3x3 13 overflow');
        disp(width);
        disp({DetM_3x3_n_13_int});
    end
	%% Проверка выходной разрядности сумматора
    if (DetM_3x3_n_13_int_width_total > width_double-1)
		disp('Sum total width in function det3x3 13 higher than 64');
        disp(width);
        disp({DetM_3x3_n_13_int});
    end
    %%
[DetM_3x3_n_14_int_sum1, DetM_3x3_n_14_int_sum1_overflow, DetM_3x3_n_14_int_sum1_abs, DetM_3x3_n_14_int_sum1_width_total] = ...
    adder(Mult_DetM_2x2_n_8_int(2),  -Mult_DetM_2x2_n_6_int(2), int_size_double, width_double);
[DetM_3x3_n_14_int, DetM_3x3_n_14_int_overflow, DetM_3x3_n_14_int_abs, DetM_3x3_n_14_int_width_total] = ...
    adder(DetM_3x3_n_14_int_sum1,  Mult_DetM_2x2_n_5_int(2), int_size_double, width_double);
    %% Проверка переполнения сумматора
    if (DetM_3x3_n_14_int_overflow == 1)
		disp('Sum det3x3 14 overflow');
        disp(width);
        disp({DetM_3x3_n_14_int});
    end
	%% Проверка выходной разрядности сумматора
    if (DetM_3x3_n_14_int_width_total > width_double-1)
		disp('Sum total width in function det3x3 12 higher than 64');
        disp(width);
        disp({DetM_3x3_n_14_int});
    end
    %%
[DetM_3x3_n_22_int_sum1, DetM_3x3_n_22_int_sum1_overflow, DetM_3x3_n_22_int_sum1_abs, DetM_3x3_n_22_int_sum1_width_total] = ...
    adder(Mult_DetM_2x2_n_10_int(1),  -Mult_DetM_2x2_n_4_int(3), int_size_double, width_double);
[DetM_3x3_n_22_int, DetM_3x3_n_22_int_overflow, DetM_3x3_n_22_int_abs, DetM_3x3_n_22_int_width_total] = ...
    adder(DetM_3x3_n_22_int_sum1,  Mult_DetM_2x2_n_3_int(3), int_size_double, width_double);
    %% Проверка переполнения сумматора
    if (DetM_3x3_n_22_int_overflow == 1)
		disp('Sum det3x3 22 overflow');
        disp(width);
        disp({DetM_3x3_n_22_int});
    end
	%% Проверка выходной разрядности сумматора
    if (DetM_3x3_n_22_int_width_total > width_double-1)
		disp('Sum total width in function det3x3 22 higher than 64');
        disp(width);
        disp({DetM_3x3_n_22_int});
    end
    %%
[DetM_3x3_n_23_int_sum1, DetM_3x3_n_23_int_sum1_overflow, DetM_3x3_n_23_int_sum1_abs, DetM_3x3_n_23_int_sum1_width_total] = ...
    adder(Mult_DetM_2x2_n_9_int(1),  -Mult_DetM_2x2_n_4_int(2), int_size_double, width_double);
[DetM_3x3_n_23_int, DetM_3x3_n_23_int_overflow, DetM_3x3_n_23_int_abs, DetM_3x3_n_23_int_width_total] = ...
    adder(DetM_3x3_n_23_int_sum1,  Mult_DetM_2x2_n_2_int(3), int_size_double, width_double);
    %% Проверка переполнения сумматора
    if (DetM_3x3_n_23_int_overflow == 1)
		disp('Sum det3x3 23 overflow');
        disp(width);
        disp({DetM_3x3_n_23_int});
    end
	%% Проверка выходной разрядности сумматора
    if (DetM_3x3_n_23_int_width_total > width_double)
		disp('Sum total width in function det3x3 23 higher than 64');
        disp(width);
        disp({DetM_3x3_n_23_int});
    end
    %%
[DetM_3x3_n_24_int_sum1, DetM_3x3_n_24_int_sum1_overflow, DetM_3x3_n_24_int_sum1_abs, DetM_3x3_n_24_int_sum1_width_total] = ...
    adder(Mult_DetM_2x2_n_8_int(1),  -Mult_DetM_2x2_n_3_int(2), int_size_double, width_double);
[DetM_3x3_n_24_int, DetM_3x3_n_24_int_overflow, DetM_3x3_n_24_int_abs, DetM_3x3_n_24_int_width_total] = ...
    adder(DetM_3x3_n_24_int_sum1,  Mult_DetM_2x2_n_2_int(2), int_size_double, width_double);
    %% Проверка переполнения сумматора
    if (DetM_3x3_n_24_int_overflow == 1)
		disp('Sum det3x3 24 overflow');
        disp(width);
        disp({DetM_3x3_n_24_int});
    end
	%% Проверка выходной разрядности сумматора
    if (DetM_3x3_n_24_int_width_total > width_double)
		disp('Sum total width in function det3x3 24 higher than 64');
        disp(width);
        disp({DetM_3x3_n_24_int});
    end
    %%
[DetM_3x3_n_33_int_sum1, DetM_3x3_n_33_int_sum1_overflow, DetM_3x3_n_33_int_sum1_abs, DetM_3x3_n_33_int_sum1_width_total] = ...
    adder(Mult_DetM_2x2_n_7_int(1),  -Mult_DetM_2x2_n_4_int(1), int_size_double, width_double);
[DetM_3x3_n_33_int, DetM_3x3_n_33_int_overflow, DetM_3x3_n_33_int_abs, DetM_3x3_n_33_int_width_total] = ...
    adder(DetM_3x3_n_33_int_sum1,  Mult_DetM_2x2_n_1_int(3), int_size_double, width_double);
    %% Проверка переполнения сумматора
    if (DetM_3x3_n_33_int_overflow == 1)
		disp('Sum det3x3 33 overflow');
        disp(width);
        disp({DetM_3x3_n_33_int});
    end
	%% Проверка выходной разрядности сумматора
    if (DetM_3x3_n_33_int_width_total > width_double)
		disp('Sum total width in function det3x3 33 higher than 64');
        disp(width);
        disp({DetM_3x3_n_33_int});
    end
    %%
[DetM_3x3_n_34_int_sum1, DetM_3x3_n_34_int_sum1_overflow, DetM_3x3_n_34_int_sum1_abs, DetM_3x3_n_34_int_sum1_width_total] = ...
    adder(Mult_DetM_2x2_n_6_int(1),  -Mult_DetM_2x2_n_3_int(1), int_size_double, width_double);
[DetM_3x3_n_34_int, DetM_3x3_n_34_int_overflow, DetM_3x3_n_34_int_abs, DetM_3x3_n_34_int_width_total] = ...
    adder(DetM_3x3_n_34_int_sum1,  Mult_DetM_2x2_n_1_int(2), int_size_double, width_double);
    %% Проверка переполнения сумматора
    if (DetM_3x3_n_34_int_overflow == 1)
		disp('Sum det3x3 34 overflow');
        disp(width);
        disp({DetM_3x3_n_34_int});
    end
	%% Проверка выходной разрядности сумматора
    if (DetM_3x3_n_34_int_width_total > width_double)
		disp('Sum total width in function det3x3 34 higher than 64');
        disp(width);
        disp({DetM_3x3_n_34_int});
    end
    %%
[DetM_3x3_n_44_int_sum1, DetM_3x3_n_44_int_sum1_overflow, DetM_3x3_n_44_int_sum1_abs, DetM_3x3_n_44_int_sum1_width_total] = ...
    adder(Mult_DetM_2x2_n_5_int(1),  -Mult_DetM_2x2_n_2_int(1), int_size_double, width_double);
[DetM_3x3_n_44_int, DetM_3x3_n_44_int_overflow, DetM_3x3_n_44_int_abs, DetM_3x3_n_44_int_width_total] = ...
    adder(DetM_3x3_n_44_int_sum1,  Mult_DetM_2x2_n_1_int(1), int_size_double, width_double);
    %% Проверка переполнения сумматора
    if (DetM_3x3_n_44_int_overflow == 1)
		disp('Sum det3x3 44 overflow');
        disp(width);
        disp({DetM_3x3_n_44_int});
    end
	%% Проверка выходной разрядности сумматора
    if (DetM_3x3_n_44_int_width_total > width_double)
		disp('Sum total width in function det3x3 44 higher than 64');
        disp(width);
        disp({DetM_3x3_n_44_int});
    end
    %%
%%
var_mult_a_det11 = a(2,2) * DetM_3x3_n_11; 
var_mult_a_det12 = a(2,3) * DetM_3x3_n_12;
var_mult_a_det13 = a(2,4) * DetM_3x3_n_13;
var_mult_a_det14 = a(2,5) * DetM_3x3_n_14;

DetM_4x4_n_1 = var_mult_a_det11 - var_mult_a_det12  + var_mult_a_det13 - var_mult_a_det14;

var_mult_a_det21 = a(2,1) * DetM_3x3_n_11;
var_mult_a_det22 = a(2,3) * DetM_3x3_n_22;
var_mult_a_det23 = a(2,4) * DetM_3x3_n_23;
var_mult_a_det24 = a(2,5) * DetM_3x3_n_24;

DetM_4x4_n_2 = var_mult_a_det21 - var_mult_a_det22 + var_mult_a_det23 - var_mult_a_det24;

var_mult_a_det31 = a(2,1) * DetM_3x3_n_12;
var_mult_a_det32 = a(2,2) * DetM_3x3_n_22;
var_mult_a_det33 = a(2,4) * DetM_3x3_n_33;
var_mult_a_det34 = a(2,5) * DetM_3x3_n_34;

DetM_4x4_n_3 = var_mult_a_det31 - var_mult_a_det32 + var_mult_a_det33 - var_mult_a_det34;

var_mult_a_det41 = a(2,1) * DetM_3x3_n_13;
var_mult_a_det42 = a(2,2) * DetM_3x3_n_23;
var_mult_a_det43 = a(2,3) * DetM_3x3_n_33;
var_mult_a_det44 = a(2,5) * DetM_3x3_n_44;

DetM_4x4_n_4 = var_mult_a_det41 - var_mult_a_det42 + var_mult_a_det43 - var_mult_a_det44;

var_mult_a_det51 = a(2,1) * DetM_3x3_n_14;
var_mult_a_det52 = a(2,2) * DetM_3x3_n_24;
var_mult_a_det53 = a(2,3) * DetM_3x3_n_34;
var_mult_a_det54 = a(2,4) * DetM_3x3_n_44;

DetM_4x4_n_5 = var_mult_a_det51 - var_mult_a_det52 + var_mult_a_det53 - var_mult_a_det54;

%%
[var_mult_a_det11_int, var_mult_a_det11_int_overflow, var_mult_a_det11_abs, var_mult_a_det11_width_total] = mult(double(a_int(2,2)), DetM_3x3_n_11_int, int_size_double, width_double);
[var_mult_a_det12_int, var_mult_a_det12_int_overflow, var_mult_a_det12_abs, var_mult_a_det12_width_total] = mult(double(a_int(2,3)), DetM_3x3_n_12_int, int_size_double, width_double);
[var_mult_a_det13_int, var_mult_a_det13_int_overflow, var_mult_a_det13_abs, var_mult_a_det13_width_total] = mult(double(a_int(2,4)), DetM_3x3_n_13_int, int_size_double, width_double);
[var_mult_a_det14_int, var_mult_a_det14_int_overflow, var_mult_a_det14_abs, var_mult_a_det14_width_total] = mult(double(a_int(2,5)), DetM_3x3_n_14_int, int_size_double, width_double);

[var_mult_a_det21_int, var_mult_a_det21_int_overflow, var_mult_a_det21_abs, var_mult_a_det21_width_total] = mult(double(a_int(2,1)), DetM_3x3_n_11_int, int_size_double, width_double);
[var_mult_a_det22_int, var_mult_a_det22_int_overflow, var_mult_a_det22_abs, var_mult_a_det22_width_total] = mult(double(a_int(2,3)), DetM_3x3_n_22_int, int_size_double, width_double);
[var_mult_a_det23_int, var_mult_a_det23_int_overflow, var_mult_a_det23_abs, var_mult_a_det23_width_total] = mult(double(a_int(2,4)), DetM_3x3_n_23_int, int_size_double, width_double);
[var_mult_a_det24_int, var_mult_a_det24_int_overflow, var_mult_a_det24_abs, var_mult_a_det24_width_total] = mult(double(a_int(2,5)), DetM_3x3_n_24_int, int_size_double, width_double);

[var_mult_a_det31_int, var_mult_a_det31_int_overflow, var_mult_a_det31_abs, var_mult_a_det31_width_total] = mult(double(a_int(2,1)), DetM_3x3_n_12_int, int_size_double, width_double);
[var_mult_a_det32_int, var_mult_a_det32_int_overflow, var_mult_a_det32_abs, var_mult_a_det32_width_total] = mult(double(a_int(2,2)), DetM_3x3_n_22_int, int_size_double, width_double);
[var_mult_a_det33_int, var_mult_a_det33_int_overflow, var_mult_a_det33_abs, var_mult_a_det33_width_total] = mult(double(a_int(2,4)), DetM_3x3_n_33_int, int_size_double, width_double);
[var_mult_a_det34_int, var_mult_a_det34_int_overflow, var_mult_a_det34_abs, var_mult_a_det34_width_total] = mult(double(a_int(2,5)), DetM_3x3_n_34_int, int_size_double, width_double);

[var_mult_a_det41_int, var_mult_a_det41_int_overflow, var_mult_a_det41_abs, var_mult_a_det41_width_total] = mult(double(a_int(2,1)), DetM_3x3_n_13_int, int_size_double, width_double);
[var_mult_a_det42_int, var_mult_a_det42_int_overflow, var_mult_a_det42_abs, var_mult_a_det42_width_total] = mult(double(a_int(2,2)), DetM_3x3_n_23_int, int_size_double, width_double);
[var_mult_a_det43_int, var_mult_a_det43_int_overflow, var_mult_a_det43_abs, var_mult_a_det43_width_total] = mult(double(a_int(2,3)), DetM_3x3_n_33_int, int_size_double, width_double);
[var_mult_a_det44_int, var_mult_a_det44_int_overflow, var_mult_a_det44_abs, var_mult_a_det44_width_total] = mult(double(a_int(2,5)), DetM_3x3_n_44_int, int_size_double, width_double);

[var_mult_a_det51_int, var_mult_a_det51_int_overflow, var_mult_a_det51_abs, var_mult_a_det51_width_total] = mult(double(a_int(2,1)), DetM_3x3_n_14_int, int_size_double, width_double);
[var_mult_a_det52_int, var_mult_a_det52_int_overflow, var_mult_a_det52_abs, var_mult_a_det52_width_total] = mult(double(a_int(2,2)), DetM_3x3_n_24_int, int_size_double, width_double);
[var_mult_a_det53_int, var_mult_a_det53_int_overflow, var_mult_a_det53_abs, var_mult_a_det53_width_total] = mult(double(a_int(2,3)), DetM_3x3_n_34_int, int_size_double, width_double);
[var_mult_a_det54_int, var_mult_a_det54_int_overflow, var_mult_a_det54_abs, var_mult_a_det54_width_total] = mult(double(a_int(2,4)), DetM_3x3_n_44_int, int_size_double, width_double);

%%
[DetM_4x4_n_1_int_sum1, DetM_4x4_n_1_int_sum1_overflow, DetM_4x4_n_1_int_sum1_abs, DetM_4x4_n_1_int_sum1_width_total] = adder(var_mult_a_det11_int,  -var_mult_a_det12_int,  int_size_double, width_double);
[DetM_4x4_n_1_int_sum2, DetM_4x4_n_1_int_sum2_overflow, DetM_4x4_n_1_int_sum2_abs, DetM_4x4_n_1_int_sum2_width_total] = adder(var_mult_a_det13_int,  -var_mult_a_det14_int,  int_size_double, width_double);
[DetM_4x4_n_1_int, 		DetM_4x4_n_1_int_overflow, 		DetM_4x4_n_1_int_abs, 		DetM_4x4_n_1_int_width_total] =		adder(DetM_4x4_n_1_int_sum1,  DetM_4x4_n_1_int_sum2, int_size_double, width_double);

[DetM_4x4_n_2_int_sum1, DetM_4x4_n_2_int_sum1_overflow, DetM_4x4_n_2_int_sum1_abs, DetM_4x4_n_2_int_sum1_width_total] = adder(var_mult_a_det21_int,  -var_mult_a_det22_int,  int_size_double, width_double);
[DetM_4x4_n_2_int_sum2, DetM_4x4_n_2_int_sum2_overflow, DetM_4x4_n_2_int_sum2_abs, DetM_4x4_n_2_int_sum2_width_total] = adder(var_mult_a_det23_int,  -var_mult_a_det24_int,  int_size_double, width_double);
[DetM_4x4_n_2_int, 		DetM_4x4_n_2_int_overflow,		DetM_4x4_n_1_int_abs, 		DetM_4x4_n_1_int_width_total] =		adder(DetM_4x4_n_2_int_sum1,  DetM_4x4_n_2_int_sum2, int_size_double, width_double);

[DetM_4x4_n_3_int_sum1, DetM_4x4_n_3_int_sum1_overflow, DetM_4x4_n_3_int_sum1_abs, DetM_4x4_n_3_int_sum1_width_total] = adder(var_mult_a_det31_int,  -var_mult_a_det32_int,  int_size_double, width_double);
[DetM_4x4_n_3_int_sum2, DetM_4x4_n_3_int_sum2_overflow, DetM_4x4_n_3_int_sum2_abs, DetM_4x4_n_3_int_sum2_width_total] = adder(var_mult_a_det33_int,  -var_mult_a_det34_int,  int_size_double, width_double);
[DetM_4x4_n_3_int, 		DetM_4x4_n_3_int_overflow,		DetM_4x4_n_3_int_abs, 		DetM_4x4_n_3_int_width_total] = 	adder(DetM_4x4_n_3_int_sum1,  DetM_4x4_n_3_int_sum2, int_size_double, width_double);

[DetM_4x4_n_4_int_sum1, DetM_4x4_n_4_int_sum1_overflow, DetM_4x4_n_4_int_sum1_abs, DetM_4x4_n_4_int_sum1_width_total] = adder(var_mult_a_det41_int,  -var_mult_a_det42_int,  int_size_double, width_double);
[DetM_4x4_n_4_int_sum2, DetM_4x4_n_4_int_sum2_overflow, DetM_4x4_n_4_int_sum2_abs, DetM_4x4_n_4_int_sum2_width_total] = adder(var_mult_a_det43_int,  -var_mult_a_det44_int,  int_size_double, width_double);
[DetM_4x4_n_4_int, 		DetM_4x4_n_4_int_overflow,		DetM_4x4_n_4_int_abs, 		DetM_4x4_n_4_int_width_total] = 	adder(DetM_4x4_n_4_int_sum1,  DetM_4x4_n_4_int_sum2, int_size_double, width_double);

[DetM_4x4_n_5_int_sum1, DetM_4x4_n_5_int_sum1_overflow, DetM_4x4_n_5_int_sum1_abs, DetM_4x4_n_5_int_sum1_width_total] = adder(var_mult_a_det51_int,  -var_mult_a_det52_int,  int_size_double, width_double);
[DetM_4x4_n_5_int_sum2, DetM_4x4_n_5_int_sum2_overflow, DetM_4x4_n_5_int_sum2_abs, DetM_4x4_n_5_int_sum2_width_total] = adder(var_mult_a_det53_int,  -var_mult_a_det54_int,  int_size_double, width_double);
[DetM_4x4_n_5_int, 		DetM_4x4_n_5_int_overflow,		DetM_4x4_n_5_int_abs, 		DetM_4x4_n_5_int_width_total] = 	adder(DetM_4x4_n_5_int_sum1,  DetM_4x4_n_5_int_sum2, int_size_double, width_double);

%%
var1 = a(1,1) * DetM_4x4_n_1;
var2 = a(1,2) * DetM_4x4_n_2;
var3 = a(1,3) * DetM_4x4_n_3;
var4 = a(1,4) * DetM_4x4_n_4;
var5 = a(1,5) * DetM_4x4_n_5;
DetM_5x5 = var1 - var2 + var3 - var4 + var5;

%%
[var1_int, var1_int_overflow, var1_int_abs, var1_int_mult_width_total] = mult(double(a_int(1,1)), DetM_4x4_n_1_int, int_size_double, width_double);
[var2_int, var2_int_overflow, var2_int_abs, var2_int_mult_width_total] = mult(double(a_int(1,2)), DetM_4x4_n_2_int, int_size_double, width_double);
[var3_int, var3_int_overflow, var3_int_abs, var3_int_mult_width_total] = mult(double(a_int(1,3)), DetM_4x4_n_3_int, int_size_double, width_double);
[var4_int, var4_int_overflow, var4_int_abs, var4_int_mult_width_total] = mult(double(a_int(1,4)), DetM_4x4_n_4_int, int_size_double, width_double);
[var5_int, var5_int_overflow, var5_int_abs, var5_int_mult_width_total] = mult(double(a_int(1,5)), DetM_4x4_n_5_int, int_size_double, width_double);

[DetM_5x5_int_sum1, DetM_5x5_int_sum1_overflow, DetM_5x5_int_sum1_abs, DetM_5x5_int_sum1_width_total] = adder(var1_int,  -var2_int, int_size_double, width_double);
[DetM_5x5_int_sum2, DetM_5x5_int_sum2_overflow, DetM_5x5_int_sum2_abs, DetM_5x5_int_sum2_width_total] = adder(var3_int,  -var4_int, int_size_double, width_double);
[DetM_5x5_int_sum3, DetM_5x5_int_sum3_overflow, DetM_5x5_int_sum3_abs, DetM_5x5_int_sum3_width_total] = adder(DetM_5x5_int_sum1,  DetM_5x5_int_sum2, 	int_size_double, width_double);
[DetM_5x5_int, DetM_5x5_int_overflow, DetM_5x5_int_abs, DetM_5x5_int_width_total] = 					adder(DetM_5x5_int_sum3,  var5_int, 			int_size_double, width_double);

end

function [det_out, det_out_int, mult1_int_abs, mult2_int_abs, det_out_int_abs, ...
    mult1_overflow, mult2_overflow, width_total_mult1, width_total_mult2, sum_overflow, width_total_sum] = det_2x2(a, a_int, int_size, width)

    det_out = a(1,1) * a(2,2) - a(2,1) * a(1,2);
	%%
    [mult1_int, mult1_overflow, mult1_int_abs, width_total_mult1] = mult(a_int(1,1), a_int(2,2), int_size, width);
    [mult2_int, mult2_overflow, mult2_int_abs, width_total_mult2] = mult(a_int(2,1), a_int(1,2), int_size, width);

	[det_out_int, sum_overflow, det_out_int_abs, width_total_sum] = adder(mult1_int,  -mult2_int, int_size, width);

end