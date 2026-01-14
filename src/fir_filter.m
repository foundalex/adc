function [y, filter_max_width_out]  = fir_filter(b, x, N, max_width_table, width, sim_options)

    max_width = table2array(max_width_table(:,2:end));

    filter_max_width_out = struct;

    buffer = cast(zeros(1,length(b)),sim_options.int_size);

	mult_n = cast(zeros(N,length(x)),sim_options.int_size);
    mult_overflow = int8(zeros(N,length(x)));
    mult_abs = cast(zeros(N,length(x)),sim_options.int_size);
    width_total_mult = int8(zeros(N,length(x)));
    filter_max_width_out.width_total_mult_max = int8(zeros(N,1));

    sum = cast(zeros(N-1,length(x)),sim_options.int_size);
	sum_overflow = int8(zeros(N-1,length(x)));
    sum_abs = cast(zeros(N-1,length(x)),sim_options.int_size);
    width_total_sum = int8(zeros(N-1,length(x)));
    filter_max_width_out.width_total_sum_max = int8(zeros(N-1,1));

    filter_max_width_out.mult_max = cast(zeros(N,1),sim_options.int_size);
    filter_max_width_out.sum_max = cast(zeros(N-1,1),sim_options.int_size);

    % buffer = [init_data]
    %% Main cycle
    for n = 1:length(x)

        buffer = cast([x(n) buffer(1:end-1)], sim_options.int_size);

		for i = uint8(1:N)
			[mult_n(i,n), mult_overflow(i,n), mult_abs(i,n), width_total_mult(i,n)] = mult(b(i), buffer(i), sim_options.int_size, width);
            %% Проверка выходной разрядности умножителей
            if (mult_overflow(i,n) == 1)
                disp(['Mult overflow in ', max_width_table.Properties.VariableNames{2}, ' number ', num2str(i)]);
                disp(width);
                disp({mult_n(i,n), i, n});
                disp(max_width_table.Properties.VariableNames{2});
            end
        
            if (width_total_mult(i,n) > width+1)
                disp(['Mult width overflow in ', max_width_table.Properties.VariableNames{2}, ' number ', num2str(i)]);
                disp(width);
                disp({width_total_mult(i,n), i, n});
                disp(max_width_table.Properties.VariableNames{2});
            end
            %% Накладываем маску
            if sim_options.enable_mask == true
                c = bitmask(mult_n(i,n), sim_options.int_size, max_width(i,1));
                if c ~= mult_n(i,n)
                    disp(['Bit mask error mult in ', max_width_table.Properties.VariableNames{2}, ' number ', num2str(i)]);
                    disp(width);
                    c = bitmask(mult_n(i,n), sim_options.int_size, max_width(i,1));
                    disp({c,mult_n(i,n)});
                    disp({sim_options.freq, sim_options.SNR});
                    disp(max_width_table.Properties.VariableNames{2});
                end
                mult_n(i,n) = c;
            else
                if filter_max_width_out.mult_max(i) < mult_abs(i,n) % определяем максимальное значение на каждом умножителе
                    filter_max_width_out.mult_max(i) = mult_abs(i,n);
                end

                % записываем макс значение разрядностей умножителей
                if filter_max_width_out.width_total_mult_max(i) < width_total_mult(i,n)
                    filter_max_width_out.width_total_mult_max(i) = width_total_mult(i,n);
                end
            end
		end
		
		%% adders
		[sum(1,n), sum_overflow(1,n), sum_abs(1,n), width_total_sum(1,n)] = adder(mult_n(1,n),  mult_n(2,n), sim_options.int_size, width);

        %% Проверка выходной разрядности сумматора
            if (sum_overflow(1,n) == 1)
                disp(['Sum overflow in ', max_width_table.Properties.VariableNames{3}, ' number 1']);
                disp(width);
                disp({sim_options.SNR, sim_options.freq});
                disp({sum(1), 1, n});
            end
            if (width_total_sum(1,n) > width+1)
                disp(['Sum width overflow in ', max_width_table.Properties.VariableNames{3}, ' number 1']);
                disp(width);
                disp({sim_options.SNR, sim_options.freq});
                disp({width_total_sum(1,n), 1, n});
            end
        %% Накладываем маску
        if sim_options.enable_mask == true
            c1 = bitmask(sum(1,n), sim_options.int_size, max_width(1,2));
            if c1 ~= sum(1,n)
                disp(['Bit mask error sum in ', max_width_table.Properties.VariableNames{3}, ' number 1']);
                disp(width);
                c1 = bitmask(sum(1,n), sim_options.int_size, max_width(1,2));
                disp({1,n});
                disp({c1,sum(1,n)});
                disp({sim_options.freq, sim_options.SNR});
            end
            sum(1,n) = c1;
        else
            % записываем макс значение на сумматоре
            if filter_max_width_out.sum_max(1) < sum_abs(1,n)
                filter_max_width_out.sum_max(1) = sum_abs(1,n);
            end

            % записываем макс значение суммы разрядностей сумматоров
            if filter_max_width_out.width_total_sum_max(1) < width_total_sum(1,n)
                filter_max_width_out.width_total_sum_max(1) = width_total_sum(1,n);
            end
        end

     	for i = uint8(1:N-2)
			[sum(i+1,n), sum_overflow(i+1,n), sum_abs(i+1,n), width_total_sum(i+1,n)] = adder(sum(i,n),  mult_n(i+2,n), sim_options.int_size, width);
            %% Проверка выходной разрядности сумматора
            if (sum_overflow(i+1,n) == 1)
                disp(['Sum overflow in ', max_width_table.Properties.VariableNames{3}, ' number ', num2str(i+1)]);
                disp(width);
                disp({sim_options.SNR, sim_options.freq});
                disp({ i+1, n});
            end
            if (width_total_sum(i+1,n) > width+1)
                disp(['Sum width overflow in ', max_width_table.Properties.VariableNames{3}, ' number ', num2str(i+1)]);
                disp(width);
                disp({sim_options.SNR, sim_options.freq});
                disp({width_total_sum(i+1,n), i+1, n});
            end
            %% Наложение маски
            if sim_options.enable_mask == true
                c1 = bitmask(sum(i+1,n), sim_options.int_size, max_width(i+1,2));
                if c1 ~= sum(i+1,n)
                    disp(['Bit mask error sum in ', max_width_table.Properties.VariableNames{3}, ' number ', num2str(i+1)]);
                    disp(width);
                    disp([i+1, n]);
                    disp([c1, sum(i+1,n)]);
                    disp([sim_options.freq, sim_options.SNR]);
                end
                sum(i+1,n) = c1;
            else
                % записываем макс значение на сумматоре
                if filter_max_width_out.sum_max(i+1) < sum_abs(i+1,n)
                    filter_max_width_out.sum_max(i+1) = sum_abs(i+1,n);
                end

                % записываем макс значение суммы разрядностей сумматоров
                if filter_max_width_out.width_total_sum_max(i+1) < width_total_sum(i+1,n)
                    filter_max_width_out.width_total_sum_max(i+1) = width_total_sum(i+1,n);
                end
            end
        end
    end

    y = sum(N-1,:)';

end