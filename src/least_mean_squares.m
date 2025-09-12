function [y_array, y_array_int, error_out] = least_mean_square(adc_input, adc_input_int, yri_cut, yri_cut_int, M, N)


% adc_input = double(adc_input_int)*2^-10;
% yri_cut = yri_cut_int*2^-11; 

figure(2)
plot([yri_cut(1:200,2), double(yri_cut_int(1:200,2))*2^-11]);


 for z = 2:M
    %% блок для расчета первых N коэффициентов фильтра
    y_out(1) = 0;
    y_out_int(1) = 0;
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

    [det_x3, det_x3_int] = determinate(double(x3_int)*2^-10, double(x3_int));
    det_func = det_x3_int * det_x3_int;
    % det_x3_int*2^-42

    det_matlab = det(x3);
    aa = det_matlab * det_matlab;
  
        %%
        for i = 1:N
            x3_shift = x3;
            x3_shift(1:N,i) = yri_cut(1:N,z);
            det_x3_shift = det(x3_shift);

            bb = det_x3_shift * det_matlab;
            www1(:,i) = bb / aa;

            %%
            x3_shift_int = x3_int;
            x3_shift_int(1:N,i) = yri_cut_int(1:N,z); 

            [det_out_shift, det_out_shift_int] = determinate(double(x3_shift_int)*2^-11, double(x3_shift_int)); % fi(1,X,55)
            % det_out_shift_int*2^-42
            det_out_mult_int = det_out_shift_int * det_x3_int; % fi(1,X,105)
            www1_int(:,i) = det_out_mult_int / det_func; 

        end

        % умножаем входные слова на рассчитанные коэффициенты
        for k = 1:N
            y_out(1) = y_out(1) + www1(k)*adc_input(k,z); % (стр 5, (13))
            % y_out = w1(1)*x(j+1) + w1(2)*x(j+2) + w1(3)*x(j+3) +
            % w1(4)*x(j+4); % Behrouz Farhang-Boroujeny, Adaptive Filters
            % Theory and Applications  (стр. 414)

            y_out_int(1) = y_out_int(1) + www1_int(k)*double(adc_input_int(k,z));
        end

        %% пересчет коэффициентов с приходом каждого слова
        for j = 1:length(yri_cut(:,1))-2*N
  
            % shift to left matrix input signal. Refresh matrix input signal
            % for every new word
            for i = 1:N
                x3(i,:) = [x3(i,2:N), 0];
                x3(i,N) = adc_input(j+N-1+i,z); % (стр.6, (20))

                x3_int(i,:) = [x3_int(i,2:N), 0];
                x3_int(i,N) = adc_input_int(j+N-1+i,z); % (стр.6, (20))
            end
        %% 

            aa = det(x3) * det(x3);

            for i = 1:N
                x3_shift = x3;
                x3_shift(1:N,i) = yri_cut(j+1:N+j,z); 
                bb = det(x3_shift) * det(x3);
                www1(:,i) = round(bb ./ aa,8);
                %%
                x3_shift_int = x3_int;
                x3_shift_int(1:N,i) = yri_cut_int(j+1:N+j,z); 

                [det_x3, det_x3_int] = determinate(double(x3_int)*2^-10, double(x3_int));
                det_func = det_x3_int * det_x3_int;
                [det_out_shift, det_out_shift_int] = determinate(double(x3_shift_int)*2^-11, double(x3_shift_int));

                det_out_mult_int = det_out_shift_int * det_x3_int;
                www1_int(:,i) = det_out_mult_int / det_func; 

            end

            % estimate coeff
            w1 = lsqminnorm(x3, yri_cut(j+1:N+j,z));
                     
            % filter input signal. Mult input words on coeff
            y_out = 0;
            y_out_int = 0;
            for k = 1:N
                y_out = y_out + www1(k)*adc_input(j+k,z); % (стр 5, (13))
                y_out_int = y_out_int + www1_int(k)*double(adc_input_int(j+k,z)); % (стр 5, (13))
            end
            y_out1(j+1) = y_out;
            y_out1_int(j+1) = y_out_int;
        end

        y_array(1:j+1,z-1) =  y_out1(1:j+1).';
        error_out(1:j+1,z-1) = y_out1(1:j+1).' ./ yri_cut(1:j+1,z);

        y_array_int(1:j+1,z-1) =  y_out1_int(1:j+1).';

 end


end


% function c = matrix_mult(a,b)
% 
% ab = length(a(:,1));
% 
%         for row = 1:ab
%             for col = 1:ab
%                 sum1 = 0;
%                 for i = 1:length(b(:,1))
%                     a1 = a(row,i);
%                     b1 = b(i,col); 
%                     sum1 = sum1 + int32(a1) * int32(b1); 
%                 end
%                 c(row,col) = sum1; 
%             end
%         end
% 
% end