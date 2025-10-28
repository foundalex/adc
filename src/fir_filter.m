function [y, mult_max, sum_max, width_total_mult_max, width_total_sum_max]  = fir_filter(b, x, enable_mask, width_mult_txt, width_sum_txt, width, sim_options)

    buffer = cast(zeros(1,length(b)),sim_options.int_size);

	mult_n = cast(zeros(sim_options.N,length(x)),sim_options.int_size);
    mult_overflow = int8(zeros(sim_options.N,length(x)));
    mult_abs = cast(zeros(sim_options.N,length(x)),sim_options.int_size);
    width_total_mult = int8(zeros(sim_options.N,length(x)));
    width_total_mult_max = int8(zeros(sim_options.N,1));

    sum = cast(zeros(sim_options.N-1,length(x)),sim_options.int_size);
	sum_overflow = int8(zeros(sim_options.N-1,length(x)));
    sum_abs = cast(zeros(sim_options.N-1,length(x)),sim_options.int_size);
    width_total_sum = int8(zeros(sim_options.N-1,length(x)));
    width_total_sum_max = int8(zeros(sim_options.N-1,1));

    mult_max = cast(zeros(sim_options.N,1),sim_options.int_size);
    sum_max = cast(zeros(sim_options.N-1,1),sim_options.int_size);

    if enable_mask == true
        width_mult = readmatrix(width_mult_txt);
        width_sum = readmatrix(width_sum_txt);
    end


    %% Main cycle
    for n = 1:length(x)

        buffer = cast([x(n) buffer(1:end-1)], sim_options.int_size);

		for i = uint8(1:sim_options.N)
			[mult_n(i,n), mult_overflow(i,n), mult_abs(i,n), width_total_mult(i,n)] = mult(b(i), buffer(i), sim_options.int_size, width);
            %% Проверка выходной разрядности умножителей
            if (mult_overflow(i,n) == 1)
                disp('Mult overflow');
                disp(width);
                disp({mult_n(i,n), i, n});
            end
        
            if (width_total_mult(i,n) > width+1)
                [mult_n(i,n), mult_overflow(i,n), mult_abs(i,n), width_total_mult(i,n)] = mult(b(i), buffer(i), sim_options.int_size, width);
                disp('Mult width overflow');
                disp(width);
                disp({width_total_mult(i,n), i, n});
            end
            %% Накладываем маску
            if enable_mask == true
                c = bitmask(mult_n(i,n), sim_options.int_size, width_mult(i));
                if c ~= mult_n(i,n)
                    disp('Bit mask error mult');
                    disp(width);
                    c = bitmask(mult_n(i,n), sim_options.int_size, width_mult(i));
                    disp({c,mult_n(i,n)});
                    disp({sim_options.freq, sim_options.SNR});
                end
                mult_n(i,n) = c;
            else
                if mult_max(i) < mult_abs(i,n) % определяем максимальное значение на кажом умножителе
                    mult_max(i) = mult_abs(i,n);
                end

                % записываем макс значение разрядностей умножителей
                if width_total_mult_max(i) < width_total_mult(i,n)
                    width_total_mult_max(i) = width_total_mult(i,n);
                end
            end
		end
		
		%% adders
		[sum(1,n), sum_overflow(1,n), sum_abs(1,n), width_total_sum(1,n)] = adder(mult_n(1,n),  mult_n(2,n), sim_options.int_size, width);

        %% Проверка выходной разрядности сумматора
            if (sum_overflow(1,n) == 1)
                disp('Sum1 overflow');
                disp(width);
                disp({sim_options.SNR, sim_options.freq});
                disp({sum(1), 1, n});
            end
            if (width_total_sum(1,n) > width+1)
                disp('Sum1 width overflow');
                disp(width);
                disp({sim_options.SNR, sim_options.freq});
                disp({width_total_sum(1,n), 1, n});
            end
        %% Накладываем маску
        if enable_mask == true
            c1 = bitmask(sum(1,n), sim_options.int_size, width_sum(1));
            if c1 ~= sum(1,n)
                disp('Bit mask error sum');
                disp(width);
                c1 = bitmask(sum(1,n), sim_options.int_size, width_sum(1));
                disp({1,n});
                disp({c1,sum(1,n)});
                disp({sim_options.freq, sim_options.SNR});
            end
            sum(1,n) = c1;
        else
            % записываем макс значение на сумматоре
            if sum_max(1) < sum_abs(1,n)
                sum_max(1) = sum_abs(1,n);
            end

            % записываем макс значение суммы разрядностей сумматоров
            if width_total_sum_max(1) < width_total_sum(1,n)
                width_total_sum_max(1) = width_total_sum(1,n);
            end
        end

     	for i = uint8(1:sim_options.N-2)
			[sum(i+1,n), sum_overflow(i+1,n), sum_abs(i+1,n), width_total_sum(i+1,n)] = adder(sum(i,n),  mult_n(i+2,n), sim_options.int_size, width);
            %% Проверка выходной разрядности сумматора
            if (sum_overflow(i+1,n) == 1)
                disp('Sum overflow');
                disp(width);
                disp({sim_options.SNR, sim_options.freq});
                disp({ i+1, n});
            end
            if (width_total_sum(i+1,n) > width+1)
                disp('Sum width overflow');
                disp(width);
                disp({sim_options.SNR, sim_options.freq});
                disp({width_total_sum(i+1,n), i+1, n});
            end
            %% Наложение маски
            if enable_mask == true
                c1 = bitmask(sum(i+1,n), sim_options.int_size, width_sum(i+1));
                if c1 ~= sum(i+1,n)
                    disp('Bit mask error sum');
                    c2 = bitmask(sum(i+1,n), sim_options.int_size, width_sum(i+1));
                    disp(width);
                    disp({i+1,n});
                    disp({c1,sum(i+1,n)});
                    disp({sim_options.freq, sim_options.SNR});
                end
                sum(i+1,n) = c1;
            else
                % записываем макс значение на сумматоре
                if sum_max(i+1) < sum_abs(i+1,n)
                    sum_max(i+1) = sum_abs(i+1,n);
                end

                % записываем макс значение суммы разрядностей сумматоров
                if width_total_sum_max(i+1) < width_total_sum(i+1,n)
                    width_total_sum_max(i+1) = width_total_sum(i+1,n);
                end
            end
        end
    end

    y = sum(72,:)';

end
