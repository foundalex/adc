function efi_out = divide(a, b, width)


    a1 = double(a);
    b1 = double(b);
    c = a1/b1;
  
    efi_out = double(c*2^width);

end
