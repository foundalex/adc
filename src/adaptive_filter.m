function [data_outd, data_out_int, adaptive_filter_structure] ...
    = adaptive_filter(input_data_double, coeff_double, input_data, coeff, adaptive_max_width, sim_options)

    % Первые 25 умножителей
    y_int = cast(zeros(sim_options.Size_matrix,sim_options.Size_matrix), sim_options.type_mult_in_adaptive_filter);
    y_int_overflow = cast(zeros(sim_options.Size_matrix,sim_options.Size_matrix), sim_options.type_mult_in_adaptive_filter);
    y_int_abs = cast(zeros(sim_options.Size_matrix,sim_options.Size_matrix), sim_options.type_mult_in_adaptive_filter);
    y_int_total_width = cast(zeros(sim_options.Size_matrix,sim_options.Size_matrix), sim_options.type_mult_in_adaptive_filter);

    y_int_out = cast(zeros(sim_options.Size_matrix*sim_options.Size_matrix,1), sim_options.type_mult_in_adaptive_filter);
    %% Первые 10 сумматоров
    sum1_data_out = cast(zeros(sim_options.Size_matrix, 2), sim_options.type_add_in_adaptive_filter);
    sum1_int_overflow = cast(zeros(sim_options.Size_matrix, 2), sim_options.type_add_in_adaptive_filter);
    sum1_array_out = cast(zeros(sim_options.Size_matrix, 2), sim_options.type_add_in_adaptive_filter);
    sum1_int_width_total = cast(zeros(sim_options.Size_matrix, 2), sim_options.type_add_in_adaptive_filter);

    sum1_int_out = cast(zeros(sim_options.Size_matrix*2,1), sim_options.type_add_in_adaptive_filter);
    %% Вторые 5 сумматоров
    sum2_data_out = cast(zeros(sim_options.Size_matrix, 1), sim_options.type_add_in_adaptive_filter);
    sum2_int_overflow = cast(zeros(sim_options.Size_matrix, 1), sim_options.type_add_in_adaptive_filter);
    %% Третьи 5 сумматоров
    sum_data_out = cast(zeros(sim_options.Size_matrix, 1), sim_options.type_add_in_adaptive_filter);
    sum_int_overflow = cast(zeros(sim_options.Size_matrix, 1), sim_options.type_add_in_adaptive_filter);
    %%
    adaptive_filter_structure = struct;

    adaptive_filter_structure.mult_int_abs = cast(zeros(sim_options.Size_matrix * sim_options.Size_matrix, 1), sim_options.type_mult_in_adaptive_filter);
    adaptive_filter_structure.mult_int_total_width = cast(zeros(sim_options.Size_matrix * sim_options.Size_matrix, 1), sim_options.type_mult_in_adaptive_filter);
    adaptive_filter_structure.sum_array_out = cast(zeros(sim_options.Size_matrix*2, 3), sim_options.type_add_in_adaptive_filter);
    adaptive_filter_structure.sum_int_width_total = cast(zeros(sim_options.Size_matrix*2,3), sim_options.type_add_in_adaptive_filter);

    kk = 1;
    %% умножители
    for i = 1:sim_options.Size_matrix
        for j = 1:sim_options.Size_matrix
            [y_int(j,i), y_int_overflow(j,i), y_int_abs(j,i), y_int_total_width(j,i)] = mult(input_data(j,i), coeff(i), sim_options.type_mult_in_adaptive_filter, sim_options.width_hilbert);   
           
            %% Проверка переполнения умножителей
            if (y_int_overflow(i) == 1)
                disp(['Mult overflow adaptive filter in ', num2str(i), num2str(j)]);
                disp({y_int_abs(i), input_data(j,i), coeff(i)});
            end
            %% Проверка выходной разрядности умножителей
            if (y_int_total_width(j,i) > sim_options.width_double)
                disp(['Mult width overflow adaptive filter in ', num2str(i), num2str(j)]);
                disp({y_int_total_width(j,i), input_data(j,i), coeff(i)});
            end

            adaptive_filter_structure.mult_int_abs(kk) = y_int_abs(j,i);
            adaptive_filter_structure.mult_int_total_width(kk) = y_int_total_width(j,i);
            y_int_out(kk) = y_int(j,i);

            %% Накладываем маску на умножители
            if sim_options.enable_mask == true
                c = bitmask(y_int_out(kk), sim_options.type_mult_in_adaptive_filter, adaptive_max_width(kk,1));
                if c ~= y_int_out(kk)
                    disp(['Bit mask error mult adaptive filter in ', num2str(kk)]);
                    c = bitmask(y_int(j,i), sim_options.type_mult_in_adaptive_filter, adaptive_max_width(kk,1));
                    disp([c, y_int_out(kk)]);
                    disp([sim_options.freq, sim_options.SNR]);
                end
                % y_int(j,i) = c;
            end

            kk = kk + 1;
        end
    end

    
    kk = 1;
    %% Первые 10 сумматоров
    for i = 1:sim_options.Size_matrix
        for j = 1:2
            [sum1_data_out(i,j), sum1_int_overflow(i,j), sum1_array_out(i,j), sum1_int_width_total(i,j)] = ...
                adder(y_int(i,j+(j-1)), y_int(i,j+j), sim_options.type_add_in_adaptive_filter, sim_options.width_hilbert);

            %% Проверка переполнения сумматоров
            if (sum1_int_overflow(i,j) == 1)
                disp(['Sum1 overflow adaptive filter in ', num2str(j), num2str(i)]);
                disp({sum1_array_out(i,j), y_int(i,j+(j-1)), y_int(i,j+j), sum1_data_out(i,j)});
            end
            %% Проверка выходной разрядности сумматоров
            if (sum1_int_width_total(i,j) > sim_options.width_hilbert)
                disp(['Sum1 width overflow adaptive filter in ', num2str(j), num2str(i)]);
                disp({sum1_array_out(i,j), y_int(i,j+(j-1)), y_int(i,j+j), sum1_data_out(i,j)});
            end

            adaptive_filter_structure.sum_array_out(kk,1) = sum1_array_out(i,j);
            adaptive_filter_structure.sum_int_width_total(kk,1) = sum1_int_width_total(i,j);
            sum1_int_out(kk) = sum1_data_out(i,j);

            %% Накладываем маску на сумматоры
            if sim_options.enable_mask == true
                c = bitmask(sum1_int_out(kk), sim_options.type_add_in_adaptive_filter, adaptive_max_width(kk,2));
                if c ~= sum1_int_out(kk)
                    disp(['Bit mask error sum1 adaptive filter in ', num2str(kk)]);
                    c = bitmask(sum1_int_out(kk), sim_options.type_add_in_adaptive_filter, adaptive_max_width(kk,2));
                    disp({c, sum1_int_out(kk)});
                    disp({sim_options.freq, sim_options.SNR});
                end
                sum1_int_out(kk) = c;
            end

            kk = kk + 1;
        end
    end

    %% Вторые 5 сумматоров
    for i = 1:sim_options.Size_matrix
        [sum2_data_out(i), sum2_int_overflow(i), adaptive_filter_structure.sum_array_out(i,2), adaptive_filter_structure.sum_int_width_total(i,2)] = ...
            adder(sum1_data_out(i,1), sum1_data_out(i,2), sim_options.type_add_in_adaptive_filter, sim_options.width_hilbert);

         %% Проверка переполнения сумматоров
        if (sum2_int_overflow(i) == 1)
            disp('Sum2 overflow adaptive filter');
            disp([adaptive_filter_structure.sum_array_out(i,2), sum2_data_out(i)]);
        end
		%% Проверка выходной разрядности сумматоров
		if (adaptive_filter_structure.sum_int_width_total(i,2) > sim_options.width_hilbert)
			disp('Sum2 width overflow adaptive filter');
			disp([adaptive_filter_structure.sum_array_out(i,2), sum2_data_out(i), sum1_data_out(i,1), sum1_data_out(i,2)]);
		end

		%% Накладываем маску на сумматоры
		if sim_options.enable_mask == true
			c = bitmask(sum2_data_out(i), sim_options.type_add_in_adaptive_filter, adaptive_max_width(i,3));
			if c ~= sum2_data_out(i)
				disp(['Bit mask error sum2 adaptive filter in ', num2str(i)]);
				c = bitmask(sum2_data_out(i), sim_options.type_add_in_adaptive_filter, adaptive_max_width(i,3));
				disp([c, sum2_data_out(i)]);
				disp([sim_options.freq, sim_options.SNR]);
			end
			sum2_data_out(i) = c;
		end
    end

	%% Третьи 5 сумматоров
    for i = 1:sim_options.Size_matrix
        [sum_data_out(i), sum_int_overflow(i), adaptive_filter_structure.sum_array_out(i,3), adaptive_filter_structure.sum_int_width_total(i,3)] = ...
            adder(sum2_data_out(i), y_int(i,5), sim_options.type_add_in_adaptive_filter, sim_options.width_hilbert);

        %% Проверка переполнения сумматоров
		if (sum_int_overflow(i) == 1)
			disp('Sum overflow adaptive filter');
			disp({adaptive_filter_structure.sum_array_out(i,3), sum_data_out(i)});
		end
		%% Проверка выходной разрядности сумматоров
		if (adaptive_filter_structure.sum_int_width_total(i,3) > sim_options.width_hilbert)
			disp('Sum width overflow adaptive filter');
			disp([adaptive_filter_structure.sum_array_out(i,3), sum_data_out(i)]);
		end

		%% Накладываем маску на сумматоры
		if sim_options.enable_mask == true
			c = bitmask(sum_data_out(i), sim_options.type_add_in_adaptive_filter, adaptive_max_width(i,4));
			if c ~= sum_data_out(i)
				disp(['Bit mask error sum adaptive filter in ', num2str(i)]);
				c = bitmask(sum_data_out(i), sim_options.type_add_in_adaptive_filter, adaptive_max_width(i,4));
				disp([c, sum_data_out(i)]);
				disp([sim_options.freq, sim_options.SNR]);
			end
			sum_data_out(i) = c;
		end
    end

    data_out_int = sum_data_out;


    %% double
    data_outd = coeff_double(1).*input_data_double(:,1)+coeff_double(2).*input_data_double(:,2)+coeff_double(3).*input_data_double(:,3)+coeff_double(4).*input_data_double(:,4)+coeff_double(5).*input_data_double(:,5);

end