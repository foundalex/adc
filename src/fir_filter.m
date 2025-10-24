function [y, mult_max, sum_max]  = fir_filter(b, x, int_size, N, enable_mask, width_mult_txt, width_sum_txt, sim_options)

    buffer = cast(zeros(1,length(b)),int_size);

	mult_n = cast(zeros(sim_options.N,length(x)),int_size);
    mult_overflow = int8(zeros(sim_options.N,length(x)));

    sum = cast(zeros(sim_options.N-1,length(x)),int_size);
	sum_overflow = int8(zeros(sim_options.N-1,length(x)));

    mult_max = cast(zeros(sim_options.N,1),int_size);
    sum_max = cast(zeros(sim_options.N-1,1),int_size);

    if enable_mask == true
        width_mult = readmatrix(width_mult_txt);
        width_sum = readmatrix(width_sum_txt);
    end

    for n = 1:length(x)

        buffer = [x(n) buffer(1:end-1)];


		for i = uint8(1:sim_options.N)
			[mult_n(i,n), mult_overflow(i,n)] = mult(b(i), cast(buffer(i), int_size), N);

            if enable_mask == true
                % Накладываем маску
                c = bitmask(mult_n(i,n), width_mult(i));
                    if c ~= mult_n(i,n)
                        disp('Bit mask error mult');
                        c = bitmask(mult_n(i,n), width_mult(i));
                        disp({c,mult_n(i,n)});
                        disp({sim_options.freq, sim_options.SNR});
                    end
                mult_n(i,n) = c;
            else
                % выясняем разрядность умножителей
                if (mult_n(i,n)) < 0
                    mult_abs = mult_n(i,n) * cast(-1, int_size); % находим число по модулю
                else
	                mult_abs = mult_n(i,n);
                end

                if mult_max(i) < mult_abs % определяем максимальное значение на кажом умножителе
                    mult_max(i) = mult_abs;
                end
            end
		end
		
		%% adders
		[sum(1,n), sum_overflow(1,n)] = adder(mult_n(1,n),  mult_n(2,n), N);

        if enable_mask == true
            c1 = bitmask(sum(1,n), width_sum(1));
            if c1 ~= sum(1,n)
                disp('Bit mask error sum');
                c1 = bitmask(sum(1,n), width_sum(1));
                disp({1,n});
                disp({c1,sum(1,n)});
                disp({sim_options.freq, sim_options.SNR});
            end
            sum(1,n) = c1;
        else
            % выясняем разрядность сумматора
            if (sum(1,n)) < 0
                sum_abs = sum(1,n) * cast(-1, int_size); % находим число по модулю
            else
	            sum_abs = sum(1,n);
            end

            if sum_max(1) < sum_abs
                sum_max(1) = sum_abs;
            end
        end

     	for i = uint8(1:sim_options.N-2)
			[sum(i+1,n), sum_overflow(i+1,n)] = adder(sum(i,n),  mult_n(i+2,n), N);
            cc = sum(i+1,n);
            if enable_mask == true
                c1 = bitmask(cc, width_sum(i+1));
                if c1 ~= cc
                    disp('Bit mask error sum');
                    c2 = bitmask(cc, width_sum(i));
                    disp({i+1,n});
                    disp({c1,sum(i+1,n)});
                    disp({sim_options.freq, sim_options.SNR});
                end
                sum(i+1,n) = c1;
            else
                % выясняем разрядность сумматоров
                if (sum(i+1,n)) < 0
                    sum_abs = sum(i+1,n) * cast(-1, int_size); % находим число по модулю
                else
	                sum_abs = sum(i+1,n);
                end

                if sum_max(i+1) < sum_abs
                    sum_max(i+1) = sum_abs;
                end
            end
        end
    end

    y = sum(72,:)';

end
