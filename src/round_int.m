function output_signal_round_int = round_int(input_signal_int, num_bit_shift, type_output)

    len = length(input_signal_int);

    output_signal_round_int = cast(zeros(len,1), type_output);

    for i = 1:len
        % округляем инты
        if (bitget(input_signal_int(i), num_bit_shift) == 1)
            output_signal_round_int(i) = bitshift(input_signal_int(i), -num_bit_shift); % сдвигаем данные
            output_signal_round_int(i) = output_signal_round_int(i) + cast(1,type_output);
        else
            output_signal_round_int(i) = bitshift(input_signal_int(i), -num_bit_shift); % сдвигаем данные
        end
    end

end