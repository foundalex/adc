
function [DetM_5x5, DetM_5x5_int] = determinate(data_in, data_in_int)

a = data_in;
a_int = data_in_int;

%% form matrix 2x2
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

[DetM_2x2_n_1, DetM_2x2_n_1_int]  = det_2x2(s.a{1}, s.a_int{1});    % a41 a42
                                                                    % a51 a52

[DetM_2x2_n_2, DetM_2x2_n_2_int]  = det_2x2(s.a{2}, s.a_int{2});    % a41 a43
                                                                    % a51 a53

[DetM_2x2_n_3, DetM_2x2_n_3_int] = det_2x2(s.a{3}, s.a_int{3});     % a41 a44
                                                                    % a41 a54

[DetM_2x2_n_4, DetM_2x2_n_4_int] = det_2x2(s.a{4}, s.a_int{4});     % a41 a45
                                                                    % a51 a55

[DetM_2x2_n_5, DetM_2x2_n_5_int] = det_2x2(s.a{5}, s.a_int{5});     % a42 a43
                                                                    % a52 a53 

[DetM_2x2_n_6, DetM_2x2_n_6_int] = det_2x2(s.a{6}, s.a_int{6});     % a42 a44
                                                                    % a52 a54

[DetM_2x2_n_7, DetM_2x2_n_7_int] = det_2x2(s.a{7}, s.a_int{7});     % a42 a45
                                                                    % a52 a55

[DetM_2x2_n_8, DetM_2x2_n_8_int] = det_2x2(s.a{8}, s.a_int{8});     % a43 a44
                                                                    % a53 a54

[DetM_2x2_n_9, DetM_2x2_n_9_int] = det_2x2(s.a{9}, s.a_int{9});     % a43 a45
                                                                    % a53 a55    

[DetM_2x2_n_10, DetM_2x2_n_10_int] = det_2x2(s.a{10}, s.a_int{10}); % a44 a45
                                                                    % a54 a55  

%% variables 1
Mult_DetM_2x2_n_10_a31 = a(3,1) * DetM_2x2_n_10;
Mult_DetM_2x2_n_10_a32 = a(3,2) * DetM_2x2_n_10;
Mult_DetM_2x2_n_10_a33 = a(3,3) * DetM_2x2_n_10;

Mult_DetM_2x2_n_7_a31 = a(3,1) * DetM_2x2_n_7;
Mult_DetM_2x2_n_7_a33 = a(3,3) * DetM_2x2_n_7;
Mult_DetM_2x2_n_7_a34 = a(3,4) * DetM_2x2_n_7;

Mult_DetM_2x2_n_9_a31 = a(3,1) * DetM_2x2_n_9;
Mult_DetM_2x2_n_9_a32 = a(3,2) * DetM_2x2_n_9;
Mult_DetM_2x2_n_9_a34 = a(3,4) * DetM_2x2_n_9;

Mult_DetM_2x2_n_8_a31 = a(3,1) * DetM_2x2_n_8;
Mult_DetM_2x2_n_8_a32 = a(3,2) * DetM_2x2_n_8;
Mult_DetM_2x2_n_8_a35 = a(3,5) * DetM_2x2_n_8;

Mult_DetM_2x2_n_6_a31 = a(3,1) * DetM_2x2_n_6;
Mult_DetM_2x2_n_6_a33 = a(3,3) * DetM_2x2_n_6;
Mult_DetM_2x2_n_6_a35 = a(3,5) * DetM_2x2_n_6;

Mult_DetM_2x2_n_5_a31 = a(3,1) * DetM_2x2_n_5;
Mult_DetM_2x2_n_5_a34 = a(3,4) * DetM_2x2_n_5;
Mult_DetM_2x2_n_5_a35 = a(3,5) * DetM_2x2_n_5;

Mult_DetM_2x2_n_4_a32 = a(3,2) * DetM_2x2_n_4;
Mult_DetM_2x2_n_4_a33 = a(3,3) * DetM_2x2_n_4;
Mult_DetM_2x2_n_4_a34 = a(3,4) * DetM_2x2_n_4;

Mult_DetM_2x2_n_3_a32 = a(3,2) * DetM_2x2_n_3;
Mult_DetM_2x2_n_3_a33 = a(3,3) * DetM_2x2_n_3;
Mult_DetM_2x2_n_3_a35 = a(3,5) * DetM_2x2_n_3;

Mult_DetM_2x2_n_2_a32 = a(3,2) * DetM_2x2_n_2;
Mult_DetM_2x2_n_2_a34 = a(3,4) * DetM_2x2_n_2;
Mult_DetM_2x2_n_2_a35 = a(3,5) * DetM_2x2_n_2;

Mult_DetM_2x2_n_1_a33 = a(3,3) * DetM_2x2_n_1;
Mult_DetM_2x2_n_1_a34 = a(3,4) * DetM_2x2_n_1;
Mult_DetM_2x2_n_1_a35 = a(3,5) * DetM_2x2_n_1;
%%
nn1 = 8;
Mult_DetM_2x2_n_10_a31_int = a_int(3,1) * DetM_2x2_n_10_int; % fi(1,12,10) * fi(1,25,20) = fi(1,37,30)
Mult_DetM_2x2_n_10_a31_int = round(Mult_DetM_2x2_n_10_a31_int/nn1); % fi(1,37,30) - 3 = fi(1,37,27) 
Mult_DetM_2x2_n_10_a32_int = a_int(3,2) * DetM_2x2_n_10_int;
Mult_DetM_2x2_n_10_a32_int = round(Mult_DetM_2x2_n_10_a32_int/nn1); 
Mult_DetM_2x2_n_10_a33_int = a_int(3,3) * DetM_2x2_n_10_int;
Mult_DetM_2x2_n_10_a33_int = round(Mult_DetM_2x2_n_10_a33_int/nn1); 

Mult_DetM_2x2_n_7_a31_int = a_int(3,1) * DetM_2x2_n_7_int;
Mult_DetM_2x2_n_7_a31_int = round(Mult_DetM_2x2_n_7_a31_int/nn1); 
Mult_DetM_2x2_n_7_a33_int = a_int(3,3) * DetM_2x2_n_7_int;
Mult_DetM_2x2_n_7_a33_int = round(Mult_DetM_2x2_n_7_a33_int/nn1); 
Mult_DetM_2x2_n_7_a34_int = a_int(3,4) * DetM_2x2_n_7_int;
Mult_DetM_2x2_n_7_a34_int = round(Mult_DetM_2x2_n_7_a34_int/nn1); 
 
Mult_DetM_2x2_n_9_a31_int = a_int(3,1) * DetM_2x2_n_9_int;
Mult_DetM_2x2_n_9_a31_int = round(Mult_DetM_2x2_n_9_a31_int/nn1); 
Mult_DetM_2x2_n_9_a32_int = a_int(3,2) * DetM_2x2_n_9_int;
Mult_DetM_2x2_n_9_a32_int = round(Mult_DetM_2x2_n_9_a32_int/nn1); 
Mult_DetM_2x2_n_9_a34_int = a_int(3,4) * DetM_2x2_n_9_int;
Mult_DetM_2x2_n_9_a34_int = round(Mult_DetM_2x2_n_9_a34_int/nn1); 

Mult_DetM_2x2_n_8_a31_int = a_int(3,1) * DetM_2x2_n_8_int;
Mult_DetM_2x2_n_8_a31_int = round(Mult_DetM_2x2_n_8_a31_int/nn1); 
Mult_DetM_2x2_n_8_a32_int = a_int(3,2) * DetM_2x2_n_8_int;
Mult_DetM_2x2_n_8_a32_int = round(Mult_DetM_2x2_n_8_a32_int/nn1); 
Mult_DetM_2x2_n_8_a35_int = a_int(3,5) * DetM_2x2_n_8_int;
Mult_DetM_2x2_n_8_a35_int = round(Mult_DetM_2x2_n_8_a35_int/nn1); 

Mult_DetM_2x2_n_6_a31_int = a_int(3,1) * DetM_2x2_n_6_int;
Mult_DetM_2x2_n_6_a31_int = round(Mult_DetM_2x2_n_6_a31_int/nn1); 
Mult_DetM_2x2_n_6_a33_int = a_int(3,3) * DetM_2x2_n_6_int;
Mult_DetM_2x2_n_6_a33_int = round(Mult_DetM_2x2_n_6_a33_int/nn1); 
Mult_DetM_2x2_n_6_a35_int = a_int(3,5) * DetM_2x2_n_6_int;
Mult_DetM_2x2_n_6_a35_int = round(Mult_DetM_2x2_n_6_a35_int/nn1); 

Mult_DetM_2x2_n_5_a31_int = a_int(3,1) * DetM_2x2_n_5_int;
Mult_DetM_2x2_n_5_a31_int = round(Mult_DetM_2x2_n_5_a31_int/nn1); 
Mult_DetM_2x2_n_5_a34_int = a_int(3,4) * DetM_2x2_n_5_int;
Mult_DetM_2x2_n_5_a34_int = round(Mult_DetM_2x2_n_5_a34_int/nn1); 
Mult_DetM_2x2_n_5_a35_int = a_int(3,5) * DetM_2x2_n_5_int;
Mult_DetM_2x2_n_5_a35_int = round(Mult_DetM_2x2_n_5_a35_int/nn1); 

Mult_DetM_2x2_n_4_a32_int = a_int(3,2) * DetM_2x2_n_4_int;
Mult_DetM_2x2_n_4_a32_int = round(Mult_DetM_2x2_n_4_a32_int/nn1); 
Mult_DetM_2x2_n_4_a33_int = a_int(3,3) * DetM_2x2_n_4_int;
Mult_DetM_2x2_n_4_a33_int = round(Mult_DetM_2x2_n_4_a33_int/nn1); 
Mult_DetM_2x2_n_4_a34_int = a_int(3,4) * DetM_2x2_n_4_int;
Mult_DetM_2x2_n_4_a34_int = round(Mult_DetM_2x2_n_4_a34_int/nn1); 

Mult_DetM_2x2_n_3_a32_int = a_int(3,2) * DetM_2x2_n_3_int;
Mult_DetM_2x2_n_3_a32_int = round(Mult_DetM_2x2_n_3_a32_int/nn1); 
Mult_DetM_2x2_n_3_a33_int = a_int(3,3) * DetM_2x2_n_3_int;
Mult_DetM_2x2_n_3_a33_int = round(Mult_DetM_2x2_n_3_a33_int/nn1); 
Mult_DetM_2x2_n_3_a35_int = a_int(3,5) * DetM_2x2_n_3_int;
Mult_DetM_2x2_n_3_a35_int = round(Mult_DetM_2x2_n_3_a35_int/nn1); 

Mult_DetM_2x2_n_2_a32_int = a_int(3,2) * DetM_2x2_n_2_int;
Mult_DetM_2x2_n_2_a32_int = round(Mult_DetM_2x2_n_2_a32_int/nn1); 
Mult_DetM_2x2_n_2_a34_int = a_int(3,4) * DetM_2x2_n_2_int;
Mult_DetM_2x2_n_2_a34_int = round(Mult_DetM_2x2_n_2_a34_int/nn1); 
Mult_DetM_2x2_n_2_a35_int = a_int(3,5) * DetM_2x2_n_2_int;
Mult_DetM_2x2_n_2_a35_int = round(Mult_DetM_2x2_n_2_a35_int/nn1); 

Mult_DetM_2x2_n_1_a33_int = a_int(3,3) * DetM_2x2_n_1_int;
Mult_DetM_2x2_n_1_a33_int = round(Mult_DetM_2x2_n_1_a33_int/nn1); 
Mult_DetM_2x2_n_1_a34_int = a_int(3,4) * DetM_2x2_n_1_int;
Mult_DetM_2x2_n_1_a34_int = round(Mult_DetM_2x2_n_1_a34_int/nn1); 
Mult_DetM_2x2_n_1_a35_int = a_int(3,5) * DetM_2x2_n_1_int;
Mult_DetM_2x2_n_1_a35_int = round(Mult_DetM_2x2_n_1_a35_int/nn1);

% mult = 5 * 10 = 50;
% add = 10;
%%
% 1 4x4
DetM_3x3_n_11 = Mult_DetM_2x2_n_10_a33 - Mult_DetM_2x2_n_9_a34 + Mult_DetM_2x2_n_8_a35;
DetM_3x3_n_12 = Mult_DetM_2x2_n_10_a32 - Mult_DetM_2x2_n_7_a34 + Mult_DetM_2x2_n_6_a35; 
DetM_3x3_n_13 = Mult_DetM_2x2_n_9_a32 - Mult_DetM_2x2_n_7_a33 + Mult_DetM_2x2_n_5_a35;
DetM_3x3_n_14 = Mult_DetM_2x2_n_8_a32 - Mult_DetM_2x2_n_6_a33 + Mult_DetM_2x2_n_5_a34;
% 2 4x4
% DetM_3x3_n_21 = Mult_DetM_2x2_n_10_a33 - Mult_DetM_2x2_n_9_a34 + Mult_DetM_2x2_n_8_a35; 
DetM_3x3_n_22 = Mult_DetM_2x2_n_10_a31 - Mult_DetM_2x2_n_4_a34 + Mult_DetM_2x2_n_3_a35; 
DetM_3x3_n_23 = Mult_DetM_2x2_n_9_a31 - Mult_DetM_2x2_n_4_a33 + Mult_DetM_2x2_n_2_a35;
DetM_3x3_n_24 = Mult_DetM_2x2_n_8_a31 - Mult_DetM_2x2_n_3_a33 + Mult_DetM_2x2_n_2_a34;
% 3 4x4
% DetM_3x3_n_31 = Mult_DetM_2x2_n_10_a32 - Mult_DetM_2x2_n_7_a34 + Mult_DetM_2x2_n_6_a35; 
% DetM_3x3_n_32 = Mult_DetM_2x2_n_10_a31 - Mult_DetM_2x2_n_4_a34 + Mult_DetM_2x2_n_3_a35; 
DetM_3x3_n_33 = Mult_DetM_2x2_n_7_a31 - Mult_DetM_2x2_n_4_a32 + Mult_DetM_2x2_n_1_a35;
DetM_3x3_n_34 = Mult_DetM_2x2_n_6_a31 - Mult_DetM_2x2_n_3_a32 + Mult_DetM_2x2_n_1_a34;
% 4 4x4
% DetM_3x3_n_41 = Mult_DetM_2x2_n_9_a32 - Mult_DetM_2x2_n_7_a33 + Mult_DetM_2x2_n_5_a35; 
% DetM_3x3_n_42 = Mult_DetM_2x2_n_9_a31 - Mult_DetM_2x2_n_4_a33 + Mult_DetM_2x2_n_2_a35; 
% DetM_3x3_n_43 = Mult_DetM_2x2_n_7_a31 - Mult_DetM_2x2_n_4_a32 + Mult_DetM_2x2_n_1_a35;
DetM_3x3_n_44 = Mult_DetM_2x2_n_5_a31 - Mult_DetM_2x2_n_2_a32 + Mult_DetM_2x2_n_1_a33;
% 5 4x4
% DetM_3x3_n_51 = Mult_DetM_2x2_n_8_a32 - Mult_DetM_2x2_n_6_a33 + Mult_DetM_2x2_n_5_a34; 
% DetM_3x3_n_52 = Mult_DetM_2x2_n_8_a31 - Mult_DetM_2x2_n_3_a33 + Mult_DetM_2x2_n_2_a34; 
% DetM_3x3_n_53 = Mult_DetM_2x2_n_6_a31 - Mult_DetM_2x2_n_3_a32 + Mult_DetM_2x2_n_1_a34;
% DetM_3x3_n_54 = Mult_DetM_2x2_n_5_a31 - Mult_DetM_2x2_n_2_a32 + Mult_DetM_2x2_n_1_a33;

% mult = 9 * 10;
% add = 5 * 10;

%%
% 1 4x4
DetM_3x3_n_11_int = Mult_DetM_2x2_n_10_a33_int - Mult_DetM_2x2_n_9_a34_int + Mult_DetM_2x2_n_8_a35_int; % fi(1,21,14) + fi(1,21,14) + fi(1,21,14) = fi(1,23,14) 
DetM_3x3_n_12_int = Mult_DetM_2x2_n_10_a32_int - Mult_DetM_2x2_n_7_a34_int + Mult_DetM_2x2_n_6_a35_int; 
DetM_3x3_n_13_int = Mult_DetM_2x2_n_9_a32_int - Mult_DetM_2x2_n_7_a33_int + Mult_DetM_2x2_n_5_a35_int;
DetM_3x3_n_14_int = Mult_DetM_2x2_n_8_a32_int - Mult_DetM_2x2_n_6_a33_int + Mult_DetM_2x2_n_5_a34_int;

DetM_3x3_n_22_int = Mult_DetM_2x2_n_10_a31_int - Mult_DetM_2x2_n_4_a34_int + Mult_DetM_2x2_n_3_a35_int; 
DetM_3x3_n_23_int = Mult_DetM_2x2_n_9_a31_int - Mult_DetM_2x2_n_4_a33_int + Mult_DetM_2x2_n_2_a35_int;
DetM_3x3_n_24_int = Mult_DetM_2x2_n_8_a31_int - Mult_DetM_2x2_n_3_a33_int + Mult_DetM_2x2_n_2_a34_int;

DetM_3x3_n_33_int = Mult_DetM_2x2_n_7_a31_int - Mult_DetM_2x2_n_4_a32_int + Mult_DetM_2x2_n_1_a35_int;
DetM_3x3_n_34_int = Mult_DetM_2x2_n_6_a31_int - Mult_DetM_2x2_n_3_a32_int + Mult_DetM_2x2_n_1_a34_int;

DetM_3x3_n_44_int = Mult_DetM_2x2_n_5_a31_int - Mult_DetM_2x2_n_2_a32_int + Mult_DetM_2x2_n_1_a33_int;
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

var_mult_a_det41 = a(2,1) * DetM_3x3_n_14;
var_mult_a_det42 = a(2,2) * DetM_3x3_n_24;
var_mult_a_det43 = a(2,3) * DetM_3x3_n_34;
var_mult_a_det44 = a(2,4) * DetM_3x3_n_44;

DetM_4x4_n_5 = var_mult_a_det41 - var_mult_a_det42 + var_mult_a_det43 - var_mult_a_det44;

% mult = 90 + 20 = 110;
% add = 50 + 15 = 66;
%%
NN = 1;

var_mult_a_det11_int = a_int(2,2) * DetM_3x3_n_11_int; % fi(1,12,10) * fi(1,X,27) = fi(1,X,37)
var_mult_a_det11_int = round(var_mult_a_det11_int/NN); % fi(1,37,19) - 1 = fi(1,X,36) 
var_mult_a_det12_int = a_int(2,3) * DetM_3x3_n_12_int;
var_mult_a_det12_int = round(var_mult_a_det12_int/NN);
var_mult_a_det13_int = a_int(2,4) * DetM_3x3_n_13_int;
var_mult_a_det13_int = round(var_mult_a_det13_int/NN);
var_mult_a_det14_int = a_int(2,5) * DetM_3x3_n_14_int;
var_mult_a_det14_int = round(var_mult_a_det14_int/NN); 

DetM_4x4_n_1_int = var_mult_a_det11_int - var_mult_a_det12_int + var_mult_a_det13_int - var_mult_a_det14_int; % fi(1,30,29) + fi(1,30,29) + fi(1,30,29) = fi(1,32,19)

var_mult_a_det21_int = a_int(2,1) * DetM_3x3_n_11_int;
var_mult_a_det21_int = round(var_mult_a_det21_int/NN);
var_mult_a_det22_int = a_int(2,3) * DetM_3x3_n_22_int;
var_mult_a_det22_int = round(var_mult_a_det22_int/NN); 
var_mult_a_det23_int = a_int(2,4) * DetM_3x3_n_23_int;
var_mult_a_det23_int = round(var_mult_a_det23_int/NN); 
var_mult_a_det24_int = a_int(2,5) * DetM_3x3_n_24_int;
var_mult_a_det24_int = round(var_mult_a_det24_int/NN);

DetM_4x4_n_2_int = var_mult_a_det21_int - var_mult_a_det22_int + var_mult_a_det23_int - var_mult_a_det24_int;

var_mult_a_det31_int = a_int(2,1) * DetM_3x3_n_12_int;
var_mult_a_det31_int = round(var_mult_a_det31_int/NN);
var_mult_a_det32_int = a_int(2,2) * DetM_3x3_n_22_int;
var_mult_a_det32_int = round(var_mult_a_det32_int/NN);
var_mult_a_det33_int = a_int(2,4) * DetM_3x3_n_33_int;
var_mult_a_det33_int = round(var_mult_a_det33_int/NN);
var_mult_a_det34_int = a_int(2,5) * DetM_3x3_n_34_int;
var_mult_a_det34_int = round(var_mult_a_det34_int/NN);

DetM_4x4_n_3_int = var_mult_a_det31_int - var_mult_a_det32_int + var_mult_a_det33_int - var_mult_a_det34_int;

var_mult_a_det41_int = a_int(2,1) * DetM_3x3_n_13_int;
var_mult_a_det41_int = round(var_mult_a_det41_int/NN);
var_mult_a_det42_int = a_int(2,2) * DetM_3x3_n_23_int;
var_mult_a_det42_int = round(var_mult_a_det42_int/NN);
var_mult_a_det43_int = a_int(2,3) * DetM_3x3_n_33_int;
var_mult_a_det43_int = round(var_mult_a_det43_int/NN);
var_mult_a_det44_int = a_int(2,5) * DetM_3x3_n_44_int;
var_mult_a_det44_int = round(var_mult_a_det44_int/NN);

DetM_4x4_n_4_int = var_mult_a_det41_int - var_mult_a_det42_int + var_mult_a_det43_int - var_mult_a_det44_int;

var_mult_a_det41_int = a_int(2,1) * DetM_3x3_n_14_int;
var_mult_a_det41_int = round(var_mult_a_det41_int/NN);
var_mult_a_det42_int = a_int(2,2) * DetM_3x3_n_24_int;
var_mult_a_det42_int = round(var_mult_a_det42_int/NN);
var_mult_a_det43_int = a_int(2,3) * DetM_3x3_n_34_int;
var_mult_a_det43_int = round(var_mult_a_det43_int/NN);
var_mult_a_det44_int = a_int(2,4) * DetM_3x3_n_44_int;
var_mult_a_det44_int = round(var_mult_a_det44_int/NN);

DetM_4x4_n_5_int = var_mult_a_det41_int - var_mult_a_det42_int + var_mult_a_det43_int - var_mult_a_det44_int;
%%
var1 = a(1,1)*DetM_4x4_n_1;
var2 = a(1,2)*DetM_4x4_n_2;
var3 = a(1,3)*DetM_4x4_n_3;
var4 = a(1,4)*DetM_4x4_n_4;
var5 = a(1,5)*DetM_4x4_n_5;
DetM_5x5 = var1 - var2 + var3 - var4 + var5;

% mult = 90 + 20 = 110 + 5 = 115;
% add = 50 + 15 = 66 + 4 + 70;
% 115 + 70 = 185
%%
r = 32;
var1_int = a_int(1,1)*DetM_4x4_n_1_int; %fi(1,X,10) * (1,X,36) = fi(1,X,46)
var1_int = round(var1_int/r);           %fi(1,X,46) - 5  = (1,X,41) 
var2_int = a_int(1,2)*DetM_4x4_n_2_int;
var2_int = round(var2_int/r);
var3_int = a_int(1,3)*DetM_4x4_n_3_int;
var3_int = round(var3_int/r);
var4_int = a_int(1,4)*DetM_4x4_n_4_int;
var4_int = round(var4_int/r);
var5_int = a_int(1,5)*DetM_4x4_n_5_int;
var5_int = round(var5_int/r);

DetM_5x5_int = var1_int - var2_int + var3_int - var4_int + var5_int; % fi(1,X,50)

%%
% a_int = data_in_int;
% 
% a1 = a(4,4)*a(5,5) - a(5,4)*a(4,5);
% a2 = a(4,3)*a(5,5) - a(5,3)*a(4,5);
% a3 = a(4,3)*a(5,4) - a(5,3)*a(4,4);
% a4 = a(4,2)*a(5,4) - a(5,2)*a(4,4);
% a5 = a(4,2)*a(5,5) - a(5,2)*a(4,5);
% a6 = a(4,2)*a(5,3) - a(5,2)*a(4,3);
% a7 = a(4,1)*a(5,5) - a(5,1)*a(4,5);
% a8 = a(4,1)*a(5,4) - a(5,1)*a(4,4);
% a9 = a(4,1)*a(5,3) - a(5,1)*a(4,3);
% a10 = a(4,1)*a(5,2) - a(5,1)*a(4,2);
% 
% aa1 = a(3,3)*a1;
% aa2 = a(3,1)*a1;
% aa3 = a(3,1)*a2;
% aa4 = a(3,2)*a1;
% aa5 = a(3,1)*a3;
% aa6 = a(3,2)*a2;
% aa7 = a(3,2)*a3;
% aa8 = a(3,1)*a5;
% aa9 = a(3,1)*a6;
% aa10 = a(3,1)*a4; 
% 
% b1 = a(3,3)*a5; c1 = a(3,5)*a3;
% b2 = a(3,3)*a4; c2 = a(3,5)*a6;
% b3 = a(3,4)*a2; c3 = a(3,4)*a6;
% b4 = a(3,4)*a7; c4 = a(3,5)*a8;
% b5 = a(3,2)*a7; c5 = a(3,5)*a9;
% b6 = a(3,3)*a7; c6 = a(3,5)*a10;
% b7 = a(3,2)*a8; c7 = a(3,3)*a10;
% b8 = a(3,2)*a9; c8 = a(3,4)*a10;
% b9 = a(3,3)*a8; c9 = a(3,4)*a9;
% b10 = a(3,4)*a5; c10 = a(3,5)*a4;
% 
% e1 = aa10 - b7 + c8;
% e2 = aa9 - b8 + c7;
% e3 = aa5 - b9 + c9;
% e4 = aa7 - b2 + c3;
% e5 = aa6 - b1 + c2;
% e6 = aa2 - b4 + c4;
% e7 = aa3 - b6 + c5;
% e8 = aa8 - b5 + c6;
% e9 = aa1 - b3 + c1;
% e10 = aa4 - b10 + c10;

% det1_4x4 = a(1,1) * (...
%            a(2,2) * e9 ...
%            -a(2,3) * e10 ...
%            +a(2,4) * e5 ...
%            -a(2,5) * e4 ...
%             );
% 
% det2_4x4 = -a(1,2) * (...
%            a(2,1) * e9 ...
%            -a(2,3) * e6 ...
%            +a(2,4) * e7 ...
%            -a(2,5) * e3 ...
%             );
% 
% det3_4x4 = a(1,3) * (...
%            a(2,1) * e10 ...
%            -a(2,2) * e6 ...
%            +a(2,4) * e8 ...
%            -a(2,5) * e1 ...
%             );
% 
% det4_4x4 = -a(1,4) * (...
%            a(2,1) * e5 ...
%            -a(2,2) * e7 ...
%            +a(2,3) * e8 ...
%            -a(2,5) * e2 ...
%             );
% 
% det5_4x4 = a(1,5) * (...
%            a(2,1) * e4 ...
%            -a(2,2) * e3 ...
%            +a(2,3) * e1 ...
%            -a(2,4) * e2 ...
%             );
% 
% sum_deet = det1_4x4 + det2_4x4 + det3_4x4 + det4_4x4 + det5_4x4;

end

function [det_out, det_out_int] = det_2x2(a, a_int)
    det_out = a(1,1)*a(2,2) - a(2,1)*a(1,2); 
    det_out_int = a_int(1,1)*a_int(2,2) - a_int(2,1)*a_int(1,2); % fi(1,12,10) * fi(1,12,10) = fi(1,24,20) + fi(1,24,20) = fi(1,25,20)
    % det_out_int = round(det_out_int/512); % fi(1,25,20) - 9 = fi(1,16,11)
end
