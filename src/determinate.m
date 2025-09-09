
function [sum_deet, sum_deet_int] = determinate(data_in, data_in_int)

a = data_in;
a_int = data_in_int; % fi(1,12,10)

a1 = a(4,4)*a(5,5) - a(5,4)*a(4,5);
a2 = a(4,3)*a(5,5) - a(5,3)*a(4,5);
a3 = a(4,3)*a(5,4) - a(5,3)*a(4,4);
a4 = a(4,2)*a(5,4) - a(5,2)*a(4,4);
a5 = a(4,2)*a(5,5) - a(5,2)*a(4,5);
a6 = a(4,2)*a(5,3) - a(5,2)*a(4,3);
a7 = a(4,1)*a(5,5) - a(5,1)*a(4,5);
a8 = a(4,1)*a(5,4) - a(5,1)*a(4,4);
a9 = a(4,1)*a(5,3) - a(5,1)*a(4,3);
a10 = a(4,1)*a(5,2) - a(5,1)*a(4,2);


a1_int = fi(a_int(4,4)*a_int(5,5) - a_int(5,4)*a_int(4,5),1,24,0); % -20
a2_int = fi(a_int(4,3)*a_int(5,5) - a_int(5,3)*a_int(4,5),1,24,0); % -20;
a3_int = fi(a_int(4,3)*a_int(5,4) - a_int(5,3)*a_int(4,4),1,24,0); % -20;
a4_int = fi(a_int(4,2)*a_int(5,4) - a_int(5,2)*a_int(4,4),1,24,0); % -20;
a5_int = fi(a_int(4,2)*a_int(5,5) - a_int(5,2)*a_int(4,5),1,24,0); % -20;
a6_int = fi(a_int(4,2)*a_int(5,3) - a_int(5,2)*a_int(4,3),1,24,0); % -20;
a7_int = fi(a_int(4,1)*a_int(5,5) - a_int(5,1)*a_int(4,5),1,24,0); % -20;
a8_int = fi(a_int(4,1)*a_int(5,4) - a_int(5,1)*a_int(4,4),1,24,0); % -20;
a9_int = fi(a_int(4,1)*a_int(5,3) - a_int(5,1)*a_int(4,3),1,24,0); % -20;
a10_int = fi(a_int(4,1)*a_int(5,2) - a_int(5,1)*a_int(4,2),1,24,0); % -20;

aa1 = a(3,3)*a1;
aa2 = a(3,1)*a1;
aa3 = a(3,1)*a2;
aa4 = a(3,2)*a1;
aa5 = a(3,1)*a3;
aa6 = a(3,2)*a2;
aa7 = a(3,2)*a3;
aa8 = a(3,1)*a5;
aa9 = a(3,1)*a6;
aa10 = a(3,1)*a4; 

aa1_int = a_int(3,3)*a1_int;
aa2_int = a_int(3,1)*a1_int;
aa3_int = a_int(3,1)*a2_int;
aa4_int = a_int(3,2)*a1_int;
aa5_int = a_int(3,1)*a3_int;
aa6_int = a_int(3,2)*a2_int;
aa7_int = a_int(3,2)*a3_int;
aa8_int = a_int(3,1)*a5_int;
aa9_int = a_int(3,1)*a6_int;
aa10_int = fi(a_int(3,1)*a4_int,1,36,30);

b1 = a(3,3)*a5; c1 = a(3,5)*a3;
b2 = a(3,3)*a4; c2 = a(3,5)*a6;
b3 = a(3,4)*a2; c3 = a(3,4)*a6;
b4 = a(3,4)*a7; c4 = a(3,5)*a8;
b5 = a(3,2)*a7; c5 = a(3,5)*a9;
b6 = a(3,3)*a7; c6 = a(3,5)*a10;
b7 = a(3,2)*a8; c7 = a(3,3)*a10;
b8 = a(3,2)*a9; c8 = a(3,4)*a10;
b9 = a(3,3)*a8; c9 = a(3,4)*a9;
b10 = a(3,4)*a5; c10 = a(3,5)*a4;

e1 = aa10 - b7 + c8;
e2 = aa9 - b8 + c7;
e3 = aa5 - b9 + c9;
e4 = aa7 - b2 + c3;
e5 = aa6 - b1 + c2;
e6 = aa2 - b4 + c4;
e7 = aa3 - b6 + c5;
e8 = aa8 - b5 + c6;
e9 = aa1 - b3 + c1;
e10 = aa4 - b10 + c10;

det1_4x4 = a(1,1) * (...
           a(2,2) * e9 ...
           -a(2,3) * e10 ...
           +a(2,4) * e5 ...
           -a(2,5) * e4 ...
            );

det2_4x4 = -a(1,2) * (...
           a(2,1) * e9 ...
           -a(2,3) * e6 ...
           +a(2,4) * e7 ...
           -a(2,5) * e3 ...
            );

det3_4x4 = a(1,3) * (...
           a(2,1) * e10 ...
           -a(2,2) * e6 ...
           +a(2,4) * e8 ...
           -a(2,5) * e1 ...
            );

det4_4x4 = -a(1,4) * (...
           a(2,1) * e5 ...
           -a(2,2) * e7 ...
           +a(2,3) * e8 ...
           -a(2,5) * e2 ...
            );

det5_4x4 = a(1,5) * (...
           a(2,1) * e4 ...
           -a(2,2) * e3 ...
           +a(2,3) * e1 ...
           -a(2,4) * e2 ...
            );

sum_deet = det1_4x4 + det2_4x4 + det3_4x4 + det4_4x4 + det5_4x4;
sum_deet_int = 0;

end