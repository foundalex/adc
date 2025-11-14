function [data_outd, data_out, y_int_abs, y_int_total_width, sum_array_out, sum_int_width_total] ...
= adaptive_filter(input_data_double, coeff_double, input_data, coeff, adaptive_mult_file, adaptive_sum_file, sim_options)

    y_int = cast(zeros(sim_options.Size_matrix,1), sim_options.type_mult_in_adaptive_filter);
    y_int_abs = cast(zeros(sim_options.Size_matrix,1), sim_options.type_mult_in_adaptive_filter);
    y_int_total_width = cast(zeros(sim_options.Size_matrix,1), sim_options.type_mult_in_adaptive_filter);

	sum_int = cast(zeros(sim_options.Size_matrix+1,1), sim_options.type_add_in_adaptive_filter);
    sum_int_abs = cast(zeros(sim_options.Size_matrix+1,1), sim_options.type_add_in_adaptive_filter);
    sum_int_width_total = cast(zeros(sim_options.Size_matrix+1,1), sim_options.type_add_in_adaptive_filter);

    if sim_options.enable_mask == true
        width_mult = readmatrix(adaptive_mult_file);
        width_sum = readmatrix(adaptive_sum_file);
    end

    sum_array_out = zeros(5,1);
    for i = 1:length(coeff)
        [y_int(i), y_int_overflow(i), y_int_abs(i), y_int_total_width(i)] = mult(input_data(i), coeff(i), sim_options.type_mult_in_adaptive_filter, sim_options.width_double);
		[sum_int(i+1), sum_int_overflow(i), sum_int_abs(i), sum_int_width_total(i)] = adder(sum_int(i), y_int(i), sim_options.type_add_in_adaptive_filter, sim_options.width_double);
        %% Проверка переполнения умножителей
        if (y_int_overflow(i) == 1)
            disp('Mult overflow adaptive filter');
            disp({y_int_abs(i), input_data(i), coeff(i)});
        end
        %% Проверка выходной разрядности умножителей
        if (y_int_total_width(i) > sim_options.width_double)
            disp('Mult width overflow adaptive filter');
            disp({y_int_total_width(i), input_data(i), coeff(i)});
        end
        %% Проверка переполнения сумматоров
        if (sum_int_overflow(i) == 1)
            disp('Sum overflow adaptive filter');
            disp({sum_int_abs(i), sum_int(i), y_int(i)});
        end
        %% Проверка выходной разрядности сумматоров
        if (y_int_total_width(i) > sim_options.width_double)
            disp('Sum width overflow adaptive filter');
            disp({y_int_total_width(i), sum_int(i), y_int(i)});
        end

        %% Накладываем маску на умножители
        if sim_options.enable_mask == true
            c = bitmask(y_int(i), sim_options.type_mult_in_adaptive_filter, width_mult(i));
             if c ~= y_int(i)
                disp('Bit mask error mult adaptivve filter');
                c = bitmask(y_int(i), sim_options.type_mult_in_adaptive_filter, width_mult(i));
                disp({c, y_int(i)});
                disp({sim_options.freq, sim_options.SNR});
            end
            y_int(i) = c;
            %% Накладываем маску на сумматоры
            c = bitmask(sum_int(i+1), sim_options.type_add_in_adaptive_filter, width_sum(i));
             if c ~= sum_int(i+1)
                disp('Bit mask error mult adaptivve filter');
                c = bitmask(sum_int(i+1), sim_options.type_add_in_adaptive_filter, width_sum(i));
                disp({c, sum_int(i+1)});
                disp({sim_options.freq, sim_options.SNR});
            end
            sum_int(i+1) = c;
        end
    end

    sum_array_out = sum_int_abs(2:6);
    data_out = sum_int(2:6);
    %%
    mult1d = input_data(1) * (coeff(1));
    mult2d = input_data(2) * (coeff(2)); 
    mult3d = input_data(3) * (coeff(3)); 
    mult4d = input_data(4) * (coeff(4)); 
    mult5d = input_data(5) * (coeff(5)); 

    data_out = mult1d + mult2d + mult3d + mult4d + mult5d;

    %% double

    mult1d = input_data_double(1) * (coeff_double(1));
    mult2d = input_data_double(2) * (coeff_double(2)); 
    mult3d = input_data_double(3) * (coeff_double(3)); 
    mult4d = input_data_double(4) * (coeff_double(4)); 
    mult5d = input_data_double(5) * (coeff_double(5)); 

    data_outd = mult1d + mult2d + mult3d + mult4d + mult5d;

end