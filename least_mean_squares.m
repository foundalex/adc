function [y_array, error_out] = least_mean_square(adc_input, yri_cut, M, N)
 for z = 2:M
        %% блок для расчета первых N коэффициентов фильтра
        y_out(1) = 0;
        % создаем матрицу входного сигнала
        for i = 1:N
            x3(i,:) = adc_input(i:N+i-1,z).'; % (стр.6, (20))
        end
        % рассчитываем первые N коэффициентов адаптивного фильтра
        % сравнивая с задержанным сигналом ADC0 (yri_cut)
        w1 = (x3'*x3) \ x3' * yri_cut(1:N,z); % (стр.6, (19))
        % w1 = lsqr(x3, adc_input_id(1:N,z));

        % умножаем входные слова на рассчитанные коэффициенты
        for k = 1:N
            y_out(1) = y_out(1) + w1(k)*adc_input(k,z); % (стр 5, (13))
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

            % estimate coeff
            w1 = ((x3'*x3) \ x3') * yri_cut(j+1:N+j,z); % (стр.6, (19))
            % w1 = lsqr(x3, adc_input_id(j+1:N+j,z));

            % filter input signal. Mult input words on coeff
            y_out = 0;
            for k = 1:N
                y_out = y_out + w1(k)*adc_input(j+k,z); % (стр 5, (13))
            end
            y_out1(j+1) = y_out;
        end

        y_array(1:j+1,z-1) =  y_out1(1:j+1).';
        error_out(1:j+1,z-1) = y_out1(1:j+1).' ./ yri_cut(1:j+1,z);

 end
end