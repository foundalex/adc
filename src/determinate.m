
function [DetM_5x5, DetM_5x5_int] = determinate(data_in, data_in_int)

a = data_in;

%% form matrix 2x2
s = struct;
e = 0;
t=4;
for n = 1:4
    for i = 1:t
        for j = 1:1
            s.a{i+e}(:,j) = a(4:end,n);
        end
        s.a{i+e}(:,2) = a(4:end,n+i);
    end
    e = e+i;
    t = t-1;
end

DetM_2x2_n_1 = det_2x2(s.a{1}); % a41 a42
                                % a51 a52
DetM_2x2_n_2 = det_2x2(s.a{2}); % a41 a43
                                % a51 a53
DetM_2x2_n_3 = det_2x2(s.a{3}); % a41 a44
                                % a41 a54
DetM_2x2_n_4 = det_2x2(s.a{4}); % a41 a45
                                % a51 a55

DetM_2x2_n_5 = det_2x2(s.a{5}); % a42 a43
                                % a52 a53 
DetM_2x2_n_6 = det_2x2(s.a{6}); % a42 a44
                                % a52 a54
DetM_2x2_n_7 = det_2x2(s.a{7}); % a42 a45
                                % a52 a55

DetM_2x2_n_8 = det_2x2(s.a{8}); % a43 a44
                                % a53 a54
DetM_2x2_n_9 = det_2x2(s.a{9}); % a43 a45
                                % a53 a55    

DetM_2x2_n_10 = det_2x2(s.a{10}); % a44 a45
                                  % a54 a55  

%%
% 1 4x4
DetM_3x3_n_11 = a(3,3) * DetM_2x2_n_10 - a(3,4) * DetM_2x2_n_9 + a(3,5) * DetM_2x2_n_8;
DetM_3x3_n_12 = a(3,2) * DetM_2x2_n_10 - a(3,4) * DetM_2x2_n_7 + a(3,5) * DetM_2x2_n_6; 
DetM_3x3_n_13 = a(3,2) * DetM_2x2_n_9 - a(3,3) * DetM_2x2_n_7 + a(3,5) * DetM_2x2_n_5;
DetM_3x3_n_14 = a(3,2) * DetM_2x2_n_8 - a(3,3) * DetM_2x2_n_6 + a(3,4) * DetM_2x2_n_5;
% 2 4x4
DetM_3x3_n_21 = a(3,3) * DetM_2x2_n_10 - a(3,4) * DetM_2x2_n_9 + a(3,5) * DetM_2x2_n_8; 
DetM_3x3_n_22 = a(3,1) * DetM_2x2_n_10 - a(3,4) * DetM_2x2_n_4 + a(3,5) * DetM_2x2_n_3; 
DetM_3x3_n_23 = a(3,1) * DetM_2x2_n_9 - a(3,3) * DetM_2x2_n_4 + a(3,5) * DetM_2x2_n_2;
DetM_3x3_n_24 = a(3,1) * DetM_2x2_n_8 - a(3,3) * DetM_2x2_n_3 + a(3,4) * DetM_2x2_n_2;
% 3 4x4
DetM_3x3_n_31 = a(3,2) * DetM_2x2_n_10 - a(3,4) * DetM_2x2_n_7 + a(3,5) * DetM_2x2_n_6; 
DetM_3x3_n_32 = a(3,1) * DetM_2x2_n_10 - a(3,4) * DetM_2x2_n_4 + a(3,5) * DetM_2x2_n_3; 
DetM_3x3_n_33 = a(3,1) * DetM_2x2_n_7 - a(3,2) * DetM_2x2_n_4 + a(3,5) * DetM_2x2_n_1;
DetM_3x3_n_34 = a(3,1) * DetM_2x2_n_6 - a(3,2) * DetM_2x2_n_3 + a(3,4) * DetM_2x2_n_1;
% 4 4x4
DetM_3x3_n_41 = a(3,2) * DetM_2x2_n_9 - a(3,3) * DetM_2x2_n_7 + a(3,5) * DetM_2x2_n_5; 
DetM_3x3_n_42 = a(3,1) * DetM_2x2_n_9 - a(3,3) * DetM_2x2_n_4 + a(3,5) * DetM_2x2_n_2; 
DetM_3x3_n_43 = a(3,1) * DetM_2x2_n_7 - a(3,2) * DetM_2x2_n_4 + a(3,5) * DetM_2x2_n_1;
DetM_3x3_n_44 = a(3,1) * DetM_2x2_n_5 - a(3,2) * DetM_2x2_n_2 + a(3,3) * DetM_2x2_n_1;
% 5 4x4
DetM_3x3_n_51 = a(3,2) * DetM_2x2_n_8 - a(3,3) * DetM_2x2_n_6 + a(3,4) * DetM_2x2_n_5; 
DetM_3x3_n_52 = a(3,1) * DetM_2x2_n_8 - a(3,3) * DetM_2x2_n_3 + a(3,4) * DetM_2x2_n_2; 
DetM_3x3_n_53 = a(3,1) * DetM_2x2_n_6 - a(3,2) * DetM_2x2_n_3 + a(3,4) * DetM_2x2_n_1;
DetM_3x3_n_54 = a(3,1) * DetM_2x2_n_5 - a(3,2) * DetM_2x2_n_2 + a(3,3) * DetM_2x2_n_1;

% mult = 9*20 = 180;
% add = 5*20 = 100
%% 
DetM_4x4_n_1 = a(2,2) * DetM_3x3_n_11 - a(2,3) * DetM_3x3_n_12 + a(2,4) * DetM_3x3_n_13 - a(2,5) * DetM_3x3_n_14;
DetM_4x4_n_2 = a(2,1) * DetM_3x3_n_21 - a(2,3) * DetM_3x3_n_22 + a(2,4) * DetM_3x3_n_23 - a(2,5) * DetM_3x3_n_24;
DetM_4x4_n_3 = a(2,1) * DetM_3x3_n_31 - a(2,2) * DetM_3x3_n_32 + a(2,4) * DetM_3x3_n_33 - a(2,5) * DetM_3x3_n_34;
DetM_4x4_n_4 = a(2,1) * DetM_3x3_n_41 - a(2,2) * DetM_3x3_n_42 + a(2,3) * DetM_3x3_n_43 - a(2,5) * DetM_3x3_n_44;
DetM_4x4_n_5 = a(2,1) * DetM_3x3_n_51 - a(2,2) * DetM_3x3_n_52 + a(2,3) * DetM_3x3_n_53 - a(2,4) * DetM_3x3_n_54;

% mult = 40*5 = 200;
% add = 5*23 = 115
%%
DetM_5x5 = a(1,1)*DetM_4x4_n_1-a(1,2)*DetM_4x4_n_2+a(1,3)*DetM_4x4_n_3-a(1,4)*DetM_4x4_n_4+a(1,5)*DetM_4x4_n_5;

% mult = 200+5 = 205;
% add = 115+4 = 119

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
DetM_5x5_int = 0;

end

function det_out = det_2x2(a)
    det_out = a(1,1)*a(2,2) - a(2,1)*a(1,2); 
end