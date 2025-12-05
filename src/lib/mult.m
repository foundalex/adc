function [y, overflow, y_abs, width_total] = mult(a, b, int_size, width)
   
    % проверка выходной разрядности
    width_a = define_of_width_int(a, int_size, width);
    width_b = define_of_width_int(b, int_size, width);

    width_total = width_a + width_b - 1;

    %% умножитель
    mmax = cast((2^(width-1))-1, int_size);
    mmin = cast(-2^(width-1), int_size);

    y = a * b;

    if (y >= mmax || y <= mmin)
        overflow = cast(1,int_size);
    else
        overflow = cast(0,int_size);
    end

    %% число по модулю
    if y < 0
        y_abs = y * cast(-1, int_size); % находим число по модулю
    else
	    y_abs = y;
    end

    % % Наложение маски на первый умножитель
    % if enable_mask == true
    %     y = bitmask(y, int_size, mask);
    % end


    % %% Проверка переполнения умножителя
    % if (overflow == 1)
	%     disp(e);
    %     % disp(sim_options.width_hilbert);
    %     % disp({i});
    % end
	% % Проверка выходной разрядности умножителя
    % if (width_total > sim_options.width_hilbert)
	%     disp(e);
    %     % disp(sim_options.width_hilbert);
    %     % disp({i});
    % end
    % 
    % % Наложение маски на первый умножитель
    % if en == true
    %     c1 = bitmask(y, sim_options.type_2x2_det, width_mult1(i));
    %     if c1 ~= y
    %         disp('Bit mask error mult 1 det 2x2');
    %         % c2 = bitmask(sum(i+1,n), sim_options.int_size, width_sum(i+1));
    %         disp(width_total);
    %         disp({c1, y});
    %         disp({sim_options.freq, sim_options.SNR});
    %     end
    %     y = c1;
    % end
end