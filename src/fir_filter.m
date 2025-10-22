function [y, min_width_mult, min_width_sum]  = fir_filter(b, x, int_size, N, enable_read, width_mult_txt, width_sum_txt)

    buffer = cast(zeros(1,length(b)),int_size);

	mult_n = cast(zeros(73,length(x)),int_size);
    mult_overflow = int8(zeros(73,length(x)));

    sum = cast(zeros(72,length(x)),int_size);
	sum_overflow = int8(zeros(72,length(x)));

    mult_max = cast(zeros(73,1),int_size);
    sum_max = cast(zeros(72,1),int_size);

    for n=1:length(x)

        buffer = [x(n) buffer(1:end-1)];

		for i = 1:73
			[mult_n(i,n), mult_overflow(i,n)] = mult(b(i), cast(buffer(i), int_size), N);

            % выясняем разрядность умножителей
            if (mult_n(i,n)) < 0
                mult_abs = mult_n(i,n) * cast(-1, int_size); % находим число по модулю
            else
	            mult_abs = mult_n(i,n);
            end

            if mult_max(i) < mult_abs
                mult_max(i) = mult_abs;
            end
		end
		
		%% adders
		[sum(1,n), sum_overflow(1,n)] = adder(mult_n(1,n),  mult_n(2,n), N);
		
     	for i = 1:71
			[sum(i+1,n), sum_overflow(i+1,n)] = adder(sum(i,n),  mult_n(i+2,n), N);

            % выясняем разрядность сумматоров
            if (sum(i,n)) < 0
                sum_abs = sum(i,n) * cast(-1, int_size); % находим число по модулю
            else
	            sum_abs = sum(i,n);
            end

            if sum_max(i) < sum_abs
                sum_max(i) = sum_abs;
            end
        end
		

    end

    y = sum(72,:)';

    if enable_read == 1
        width_mult = cast(readmatrix(width_mult_txt), int_size);
        width_sum = cast(readmatrix(width_sum_txt), int_size);
    else
        width_mult = cast(define_of_width_int(mult_max), int_size);
        width_sum = cast(define_of_width_int(sum_max), int_size);
    end
    
    %%

    % for i = 1:length(mult_n(:,1))
    %     for j = 1:length(mult_n(1,:))
    %         [c,a] = bitmask(mult_n(i,j), width_mult(i), int_size);
    %         if a ~= c
    %             disp('Bit mask error mult');
    %             disp({i,j});
    %         end
    %     end
    % end
    % 
    % for i = 1:length(sum(:,1))
    %     for j = 1:length(sum(1,:))
    %         [c1,a1] = bitmask(sum(i,j), width_sum(i), int_size);
    %         if a1 ~= c1
    %             disp('Bit mask error sum');
    %             disp({i,j});
    %         end
    % 
    %     end
    % end

    min_width_mult = 0; % cast(width_mult, int_size);
    min_width_sum = 0; %cast(width_sum, int_size);
	
end
