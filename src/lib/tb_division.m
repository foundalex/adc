clear all;
clc;
%%
% for i = 1:1000

    dat = randi([-16 16],1,2); % int

    aa = fi(dat,1,5,0);
    aa = fi(7,1,5,0);
    bb = fi(16,1,6,0);
    out = int_division(aa, bb, 5);

    double(out) * 2^-17
    dat(1)/dat(2)