function [data_outd, data_out, y_int_abs, y_int_total_width, sum_array_out, sum_int_width_total] ...
= adaptive_filter(input_data_double, coeff_double, input_data, coeff, adaptive_mult_file, adaptive_sum_file, sim_options)

    y_int = cast(zeros(sim_options.Size_matrix,sim_options.Size_matrix), sim_options.type_mult_in_adaptive_filter);
    y_int_abs = cast(zeros(sim_options.Size_matrix,sim_options.Size_matrix), sim_options.type_mult_in_adaptive_filter);
    y_int_total_width = cast(zeros(sim_options.Size_matrix,sim_options.Size_matrix), sim_options.type_mult_in_adaptive_filter);

	data_out = cast(zeros(sim_options.Size_matrix,1), sim_options.type_add_in_adaptive_filter);
    sum_array_out = cast(zeros(sim_options.Size_matrix,1), sim_options.type_add_in_adaptive_filter);
    sum_int_width_total = cast(zeros(sim_options.Size_matrix,1), sim_options.type_add_in_adaptive_filter);

    if sim_options.enable_mask == true
        width_mult = readmatrix(adaptive_mult_file);
        width_sum = readmatrix(adaptive_sum_file);
    end


    %% умножители
    for i = 1:sim_options.Size_matrix
        for j = 1:sim_options.Size_matrix
            [y_int(j,i), y_int_overflow(j,i), y_int_abs(j,i), y_int_total_width(j,i)] = mult(input_data(j,i), coeff(i), sim_options.type_mult_in_adaptive_filter, sim_options.width_hilbert);   
           
            %% Проверка переполнения умножителей
            if (y_int_overflow(i) == 1)
                disp('Mult overflow adaptive filter');
                disp({y_int_abs(i), input_data(j,i), coeff(i)});
            end
            %% Проверка выходной разрядности умножителей
            if (y_int_total_width(j,i) > sim_options.width_double)
                disp('Mult width overflow adaptive filter');
                disp({y_int_total_width(j,i), input_data(j,i), coeff(i)});
            end

            %% Накладываем маску на умножители
            if sim_options.enable_mask == true
                c = bitmask(y_int(j,i), sim_options.type_mult_in_adaptive_filter, width_mult(i));
                if c ~= y_int(j,i)
                    disp('Bit mask error mult adaptivve filter');
                    c = bitmask(y_int(j,i), sim_options.type_mult_in_adaptive_filter, width_mult(i));
                    disp({c, y_int(j,i)});
                    disp({sim_options.freq, sim_options.SNR});
                end
                y_int(j,i) = c;
            end
        end
    end

    %% сумматоры
    for i = 1:sim_options.Size_matrix
        for j = 1:sim_options.Size_matrix
            [data_out(i), sum_int_overflow(i), sum_array_out(i), sum_int_width_total(i)] = adder(data_out(i), y_int(i,j), sim_options.type_add_in_adaptive_filter, sim_options.width_hilbert);

            %% Проверка переполнения сумматоров
            if (sum_int_overflow(i) == 1)
                disp('Sum overflow adaptive filter');
                disp({sum_array_out(i), data_out(i), y_int(i,j)});
            end
            %% Проверка выходной разрядности сумматоров
            if (y_int_total_width(i) > sim_options.width_double)
                disp('Sum width overflow adaptive filter');
                disp({sum_array_out(i), data_out(i), y_int(i,j)});
            end

            %% Накладываем маску на сумматоры
            if sim_options.enable_mask == true
                c = bitmask(data_out(i), sim_options.type_add_in_adaptive_filter, width_sum(i));
                if c ~= data_out(i)
                    disp('Bit mask error mult adaptivve filter');
                    c = bitmask(data_out(i), sim_options.type_add_in_adaptive_filter, width_sum(i));
                    disp({c, data_out(i)});
                    disp({sim_options.freq, sim_options.SNR});
                end
                data_out(i) = c;
            end
        end
    end

    %% double
    data_outd = coeff_double(1).*input_data_double(:,1)+coeff_double(2).*input_data_double(:,2)+coeff_double(3).*input_data_double(:,3)+coeff_double(4).*input_data_double(:,4)+coeff_double(5).*input_data_double(:,5);

end