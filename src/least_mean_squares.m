function [y_array, error_out] = least_mean_square(adc_input, adc_input_int, yri_cut, yri_cut_int, M, N)

 for z = 2:M
    %% блок для расчета первых N коэффициентов фильтра
    y_out(1) = 0;
    % создаем матрицу входного сигнала
    for i = 1:N
        x3(i,:) = adc_input(i:N+i-1,z).'; % (стр.6, (20))
        x3_int(i,:) = adc_input_int(i:N+i-1,z).'; % fi(1,12,10)
    end
    % рассчитываем первые N коэффициентов адаптивного фильтра
    % сравнивая с задержанным сигналом ADC0 (yri_cut)
    % w1 = (x3' * x3) \ x3' * yri_cut(1:N,z); % (стр.6, (19))
    % w1 = linsolve(x3, yri_cut(1:N,z));
    % w1 = lsqr(x3, adc_input_id(1:N,z));

    w1 = lsqminnorm(x3, yri_cut(1:N,z));

    [aa1, aa_int1] = determinate(x3, double(x3_int));
    [aa2, aa_int2] = determinate(x3, double(x3_int));

    [ee, ee1] = determinate(x3*x3);

    bb = aa1*aa2;

    aa = det(x3) * det(x3);
    aa_int = aa_int1 * aa_int2; % fi (1,X,100)
    det_a_int_in_double = aa_int*2^-100;
  
        %%
        for i = 1:N
            x3_shift = x3;
            x3_shift(1:N,i) = yri_cut(1:N,z);
            bb = det(x3_shift) * det(x3);
            www1(:,i) = round(bb / aa,8);
            %%

            x3_shift_int = x3_int;
            x3_shift_int(1:N,i) = yri_cut_int(1:N,z); 
            
            [det_out, det_out_int] = determinate(x3, double(x3_int));
            [det_out_shift, det_out_shift_int] = determinate(x3_shift, double(x3_shift_int));


            det_out_mult = det_out * det_out_shift;
            det_out_mult_int = det_out_int * det_out_shift_int;
            www1_int(:,i) = det_out_mult_int / aa_int; % fi(1,X,1)

        end

        % figure(4);
        % plot([w1(:,1), double(w1_int_part1(:,1)) * 2^-20, double(w1_int_part1_int(:,1)) * 2^-7]);
        %%

        % умножаем входные слова на рассчитанные коэффициенты
        for k = 1:N
            y_out(1) = y_out(1) + www1(k)*adc_input(k,z); % (стр 5, (13))
            % y_out = w1(1)*x(j+1) + w1(2)*x(j+2) + w1(3)*x(j+3) +
            % w1(4)*x(j+4); % Behrouz Farhang-Boroujeny, Adaptive Filters
            % Theory and Applications  (стр. 414)
        end

        %% пересчет коэффициентов с приходом каждого слова
        for j = 1:length(yri_cut(:,1))-2*N
  
            % shift to left matrix input signal. Refresh matrix input signal
            % for every new word
            for i = 1:N
                x3(i,:) = [x3(i,2:N), 0];
                x3(i,N) = adc_input(j+N-1+i,z); % (стр.6, (20))
            end
        %% 

            aa = det(x3) * det(x3');

            for i = 1:N
                x3_shift = x3;
                x3_shift(1:N,i) = yri_cut(j+1:N+j,z); 
                bb = det(x3_shift) * det(x3);
                www1(:,i) = round(bb ./ aa,8);
            end

            % estimate coeff
            w1 = lsqminnorm(x3, yri_cut(j+1:N+j,z));
                     
            % filter input signal. Mult input words on coeff
            y_out = 0;
            for k = 1:N
                y_out = y_out + www1(k)*adc_input(j+k,z); % (стр 5, (13))
            end
            y_out1(j+1) = y_out;
        end

        y_array(1:j+1,z-1) =  y_out1(1:j+1).';
        error_out(1:j+1,z-1) = y_out1(1:j+1).' ./ yri_cut(1:j+1,z);

 end
end


function c = matrix_mult(a,b)

ab = length(a(:,1));

        for row = 1:ab
            for col = 1:ab
                sum1 = 0;
                for i = 1:length(b(:,1))
                    a1 = a(row,i);
                    b1 = b(i,col); 
                    sum1 = sum1 + int32(a1) * int32(b1); 
                end
                c(row,col) = sum1; 
            end
        end

end