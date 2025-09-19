function [y_array, y_array_int, error_out, error_det, error_det_lu] = least_mean_square(adc_input, adc_input_int, yri_cut, yri_cut_int, M, N)


% adc_input = double(adc_input_int)*2^-10;
% yri_cut = yri_cut_int*2^-11; 

% figure(2)
% plot([yri_cut(1:200,2), double(yri_cut_int(1:200,2))*2^-11]);

kk = 0;
tt = 0;

 for z = 2:M
    %% блок для расчета первых N коэффициентов фильтра
    y_out(1) = 0;
    % y_out_int = 0;
    % создаем матрицу входного сигнала
    for i = 1:N
        x3(i,:) = adc_input(i:N+i-1,z).'; % (стр.6, (20))
        x3_int(i,:) = adc_input_int(i:N+i-1,z).'; % fi(1,12,11)
    end
    % рассчитываем первые N коэффициентов адаптивного фильтра
    % сравнивая с задержанным сигналом ADC0 (yri_cut)
    % w1 = (x3' * x3) \ x3' * yri_cut(1:N,z); % (стр.6, (19))
    % w1 = linsolve(x3, yri_cut(1:N,z));
    % w1 = lsqr(x3, adc_input_id(1:N,z));

    % w1 = lsqminnorm(x3, yri_cut(1:N,z));
    % w1_int = lsqminnorm(double(x3_int)*2^-11, double(yri_cut_int(1:N,z))*2^-11);
   
    [det_x3, det_x3_int] = determinate(double(x3_int)*2^-11, x3_int);
    det_func = det_x3_int;


    tt = tt + 1;
    det_matlab(tt) = det(x3);

    aa = det_matlab(tt) * det_matlab(tt);

        %%
        for i = 1:N
            kk = kk + 1;

            x3_shift = x3;
            x3_shift(1:N,i) = yri_cut(1:N,z);
            
            det_x3_shift(kk) = det(x3_shift);

            bb = det_x3_shift(kk) * det_matlab(tt);
            www1(:,i) = bb / aa;

            %%
            x3_shift_int = x3_int;
            x3_shift_int(1:N,i) = yri_cut_int(1:N,z); 

            [det_out_shift, det_out_shift_int] = determinate(double(x3_shift_int)*2^-11, x3_shift_int);
            det_out_mult_int(i) = det_out_shift_int;

            error_det(kk) = abs(double(det_out_shift_int)*2^-55 / det_out_shift);
            error_det_lu(kk) = abs(det_out_shift/ det(double(x3_shift_int)*2^-11));
            %% divide
            if (det_out_mult_int(i) == 0)
                det_out_mult_int(i) = fi(1,1,70,0);
            end
            if (det_func == 0)
                det_func = fi(1,1,70,0);
            end

            % www1_int(:,i) = double(det_out_mult_int(i)) / double(det_func); 
            % T = numerictype('Signed', true,'WordLength', 70, 'FractionLength', 55);
            % www1_int(:,i) = divide(T, det_out_mult_int(i), det_func);

            www1_int(:,i) = int_division(det_out_mult_int(i), det_func, 55);
        end

            % e(i) = fi(double(www1_int(:,i))*2^-55,1,70,55);
            % zz(kk) = (double(e(i))) / (double(www1_int11(:,i)));

        % умножаем входные слова на рассчитанные коэффициенты
        for k = 1:N
            y_out(1) = y_out(1) + www1(k) * adc_input(k,z); % (стр 5, (13))
            % y_out = w1(1)*x(j+1) + w1(2)*x(j+2) + w1(3)*x(j+3) +  w1(4)*x(j+4); % Behrouz Farhang-Boroujeny, Adaptive Filters Theory and Applications  (стр. 414)
            dat_in_filt(k) = fi(adc_input_int(k,z),1,12,0);
            % y_out_int(1) = y_out_int(1) + www1_int(k) * fi(adc_input_int(k,z),1,12,0);
        end

        y_out_int = filter_transversal(dat_in_filt, www1_int);

        y_out1_int(1) = y_out_int;
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
            % estimate coeff
            % w1 = lsqminnorm(x3, yri_cut(j+1:N+j,z));

            tt = tt + 1;
            det_matlab(tt) = det(x3);
            aa = det_matlab(tt) * det_matlab(tt);

            for i = 1:N
                kk = kk + 1;

                x3_shift = x3;
                x3_shift(1:N,i) = yri_cut(j+1:N+j,z); 

                det_x3_shift(kk) = det(x3_shift);

                bb = det_x3_shift(kk) * det(x3);
                www1(:,i) = bb ./ aa;
                %%
                x3_shift_int = x3_int;
                x3_shift_int(1:N,i) = yri_cut_int(j+1:N+j,z); 

                [det_x3, det_x3_int] = determinate(double(x3_int)*2^-11, x3_int);
                det_func(i) = det_x3_int;
                [det_out_shift, det_out_shift_int] = determinate(double(x3_shift_int)*2^-11, x3_shift_int);
                det_out_mult_int(i) = det_out_shift_int;

                if (j==17 & i ==3)
                    eq = 1;
                end

                %% divide
                if (det_out_mult_int(i) == 0)
                    det_out_mult_int(i) = fi(1,1,70,0);
                end

                if (det_func(i) == 0)
                    det_func(i) = fi(1,1,70,0);
                end

                % www1_int(:,i) = double(det_out_mult_int(i)) / double(det_func(i)); 
                % T = numerictype('Signed', true,'WordLength', 70, 'FractionLength', 55);
                % www1_int(:,i) = divide(T, det_out_mult_int(i), det_func(i));

                www1_int(:,i) = int_division(det_out_mult_int(i), det_func(i), 55);

                % error_det(kk) = abs(((double(det_out_shift_int)*2^-55) / det_out_shift));
                % error_det_lu(kk) = abs(det_out_shift/ det(double(x3_shift_int)*2^-11));

                % if (isinf(www1_int(:,i)) | isnan(www1_int(:,i)))
                %     w = 1;
                % end
                % e(i) = fi(double(www1_int(:,i))*2^-55,1,70,55);
                % zz(kk) = (double(e(i))) / (double(www1_int11(:,i)));      

                y_out = 0;
                % y_out_int = 0;

                % filter input signal. Mult input words on coeff
                for k = 1:N
                    y_out = y_out + www1(k) * adc_input(j+k,z); % (стр 5, (13))

                    % y_out_int = y_out_int + www1_int(k) * double(adc_input_int(j+k,z))*2^-11;
                    dat_in_filt(k) = fi(adc_input_int(j+k,z),1,12,0);
                end

                y_out_int = filter_transversal(dat_in_filt, www1_int);

            end
            y_out1(j+1) = y_out;

            if (abs(y_out_int) > 1*10^23)
                eq = 1;
            end

            y_out1_int(j+1) = y_out_int;
        end

        y_array(1:j+1,z-1) =  y_out1(1:j+1).';
        error_out(1:j+1,z-1) = y_out1(1:j+1).' ./ yri_cut(1:j+1,z);

        y_array_int(1:j+1,z-1) = y_out1_int(1:j+1).';

 end

 min_det_matrix = min(det_matlab);
 min_shift_det_matrix = min(det_x3_shift); 

 max_det_matrix = max(det_matlab);
 max_shift_det_matrix = max(det_x3_shift);

 minimum = min([min_det_matrix min_shift_det_matrix]);
 maximum = max([max_det_matrix max_shift_det_matrix]);

end