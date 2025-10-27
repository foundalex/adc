function width = define_of_width_int(input, int_size, num)

width = int8(0);
one = int8(1);

divisor = cast(2,int_size);

%% находим модуль числа
if (input < 0)
    input_abs = input * cast(-1, int_size);
else
    input_abs = input;
end
%%

% quotient = input_abs / divisor;
% width = width + one;

for k = 1:num
    if input_abs == 1 | input_abs == 0
        % width = width + one;
        break;
    else
        input_abs = input_abs / divisor;
        width = width + one;
    end
end

width = width + one; % добавляем 1 разряд для знака

end