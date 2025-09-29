function efi_out = int_division(a, b, n, width, length_word)


    a1 = double(a);
    b1 = double(b);
    c = a1/b1;

    % T = numerictype('Signed', true,'WordLength', n, 'FractionLength', 0);
    % www = divide(T, a, b);
    % www_fix_double = double(www);
    % 
    % 
    % m = fi((2^70)-1,1,70+1,0);
    % dividend = m*a;
    % out5 = (dividend/b)*2^-55;
    % out6 = double(out5);
    % 
    % 
    % m1 = (2^(n));
    % 
    % dividend1 = (a * m1);

    % T1 = numerictype('Signed', true,'WordLength', 125, 'FractionLength', 0);
    % out = divide(T1, dividend1, b);
    % 
    % outf = fi(out,1,70,0);
    % outd = double(outf)*2^-70;


    ee = double(a)/double(b);
    e1 = ee;


    efi = fi(e1,1,length_word,width);
    efi_d = double(efi)*2^width;
    efi_out = fi(efi_d,1,length_word,0);


    out = e1*2^width;
    out1 = fi(out,1,n,0);
    out2 = double(out1)*2^-width;


   


    % enq = fi(out1,1,72,0);
    % out2 = bitshift(enq,-52);
    % out3 = fi(out2,1,20,0);
    % 
    % out3_double = double(out3)*2^-18;

    % most_significant_bit = dec2bin(bitget(out1,86:-1:71));

    % a1 = strcat(string('0b'), num2str(most_significant_bit).', string('s16'));
    % str = convertStringsToChars(a1);
    % a2 = bin2dec(str);
    % 
    % quotient1 = fi(a2,1,15,0);
    % quotient_double = double(quotient1);

    %% fractional
    % lsb = dec2bin(bitget(out1,70:-1:55));
    % a3 = strcat(string('0b'), num2str(lsb).', string('s16'));
    % str1 = convertStringsToChars(a3);
    % a4 = bin2dec(str1);

end
