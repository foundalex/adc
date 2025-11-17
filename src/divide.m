function [efi_out, overflow, c_abs, width_total] = divide(a, b, int_size_dividend, int_size_divisor, width_dividend, width_divisor, int_size_quotient, factor)

    % проверка выходной разрядности
    width_a = define_of_width_int(a, int_size_dividend, width_dividend);
    width_b = define_of_width_int(b, int_size_divisor, width_divisor);

    if (width_a > width_b)
        width_total = width_a;
    else
        width_total = width_b;
    end

    %% делитель
    mmax = cast((2^(width_divisor-1))-1, int_size_dividend);
    mmin = cast(-2^(width_divisor-1), int_size_dividend);

    a1 = double(a);
    b1 = double(b);
    c = a1 / b1;
  
    efi_out = cast(round(double(c*2^factor)),int_size_quotient);
    % efi_out = c;

    if (efi_out >= mmax || efi_out <= mmin)
        overflow = cast(1,int_size_quotient);
    else
        overflow = cast(0,int_size_quotient);
    end
    %% число по модулю
    if c < 0
        c_abs = efi_out * cast(-1, int_size_quotient); % находим число по модулю
    else
	    c_abs = efi_out;
    end
end
