function efi_out = int_division(a, b, length_word, width)


    a1 = double(a);
    b1 = double(b);
    c = a1/b1;

    % T = numerictype('Signed', true,'WordLength', n, 'FractionLength', 0);
    % www = divide(T, a, b);
    % www_fix_double = double(www);
  

    ee = double(a)/double(b);
    e1 = ee;

    efi = fi(e1,1,length_word,width);
    efi_d = double(efi)*2^width;
    efi_out = fi(efi_d,1,length_word,0);


    out = e1*2^width;
    out1 = fi(out,1,length_word,0);
    out2 = double(out1)*2^-width;

end
