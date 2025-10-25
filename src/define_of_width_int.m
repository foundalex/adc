function width = define_of_width_int(input,int_size)

width = int8(0);
one = int8(1);

divisor = cast(2,int_size);

quotient = input / divisor;
width = width + one;

for k = 1:64
    if quotient == 1
        width = width + one;
        break;
    elseif (quotient == 0) 
        break;
    else
        width = width + one;
        quotient = quotient / divisor;
    end
end

width = width + one; % добавляем 1 разряд для знака

end