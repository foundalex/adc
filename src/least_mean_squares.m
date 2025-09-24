function [y_array, y_array_int] = least_mean_square(adc_input, adc_input_int, yri_cut, yri_cut_int, M, N)

kk = 0;
tt = 0;
width = 14;
n = 70;

 for z = 2:M
    %% блок для расчета первых N коэффициентов фильтра

    % создаем матрицу входного сигнала
    for i = 1:N
        x3(i,:) = adc_input(i:N+i-1,z).'; % (стр.6, (20))
        x3_int(i,:) = adc_input_int(i:N+i-1,z).'; % fi(1,12,11)
    end

    % рассчитываем первые N коэффициентов адаптивного фильтра
    % сравнивая с задержанным сигналом ADC0 (yri_cut)
    % w1 = (x3' * x3) \ x3' * yri_cut(1:N,z); % (стр.6, (19))
    % w1 = lsqminnorm(x3, yri_cut(1:N,z));
    % w1_int = lsqminnorm(double(x3_int)*2^-11, double(yri_cut_int(1:N,z))*2^-11);
   
    [det_x3, det_x3_int] = determinate(x3, x3_int);

    tt = tt + 1;
    det_matlab(tt) = det(x3);

        %%
        for i = 1:N
            kk = kk + 1;

            x3_shift = x3;
            x3_shift(1:N,i) = yri_cut(1:N,z);
            
            det_x3_shift(kk) = det(x3_shift);

            if (det_matlab(tt) == 0)
                det_matlab(tt) = 1;
            end

            if (det_x3_shift(kk) == 0)
                det_x3_shift(kk) = 1;
            end

            www1(:,i) = det_x3_shift(kk) / det_matlab(tt);

            %% integer determinant
            x3_shift_int = x3_int;
            x3_shift_int(1:N,i) = yri_cut_int(1:N,z); 

            [det_out_shift, det_out_shift_int] = determinate(x3_shift, x3_shift_int);

            %% integer divide
            if (det_out_shift_int == 0)
                det_out_shift_int = fi(1,1,70,0);
            end
            if (det_x3_int == 0)
                det_x3_int = fi(1,1,70,0);
            end

            www1_int(:,i) = int_division(det_out_shift_int, det_x3_int, n, width);
            www1_int_double(i,:) = double(www1_int(:,i))*2^-14;
        end

        % умножаем входные слова на рассчитанные коэффициенты
        % y_out = w1(1)*x(j+1) + w1(2)*x(j+2) + w1(3)*x(j+3) +  w1(4)*x(j+4); % Behrouz Farhang-Boroujeny, Adaptive Filters Theory and Applications  (стр. 414)

        for k = 1:N
            dat_in_filt_double(k) = adc_input(k,z);
            dat_in_filt(k) = fi(adc_input_int(k,z),1,12,0);
        end

        [y_out, y_out_int] = filter_transversal(dat_in_filt_double, www1, dat_in_filt, www1_int);

        y_array(1,z-1) = y_out;
        y_array_int(1,z-1) = y_out_int;


        %% part 2

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
            tt = tt + 1;
            det_matlab(tt) = det(x3);

            if (j==17 & i ==5 & tt ==18)
                eq = 1;
            end
            [det_x3(tt), det_x3_int(tt)] = determinate(double(x3_int)*2^-11, x3_int);

            er_det_x3(tt) = det_matlab(tt) / det_x3(tt); 
            er_det_x3_1(tt) = det_matlab(tt) / (double(det_x3_int(tt))*2^-55);


            for i = 1:N
                kk = kk + 1;

                x3_shift = x3;
                x3_shift(1:N,i) = yri_cut(j+1:N+j,z); 

                det_x3_shift(kk) = det(x3_shift);

                if (det_matlab(tt) == 0)
                    det_matlab(tt) = 1;
                end

                if (det_x3_shift(kk) == 0)
                    det_x3_shift(kk) = 1;
                end

                www1(:,i) = det_x3_shift(kk) ./ det_matlab(tt);

                %% determinant
                x3_shift_int = x3_int;
                x3_shift_int(1:N,i) = yri_cut_int(j+1:N+j,z); 

                [det_out_shift, det_out_shift_int] = determinate(double(x3_shift_int)*2^-11, x3_shift_int);

                %% divide
                if (det_out_shift_int == 0)
                    det_out_shift_int = fi(1,1,70,0);
                end

                if (det_x3_int(tt) == 0)
                    det_x3_int(tt) = fi(1,1,70,0);
                end

                if (j == 401 & i == 5 & z == 3)
                    w = 1;
                end  

                www1_int(:,i) = int_division(det_out_shift_int, det_x3_int(tt), n, width);
                www1_int_double(i,:) = double(www1_int(:,i))*2^-14;

                %% filter

                % y_out = 0;

                % filter input signal. Mult input words on coeff
                for k = 1:N
                    % y_out = y_out + www1(k) * adc_input(j+k,z); % (стр 5, (13))

                    dat_in_filt_double(k) = adc_input(j+k,z);
                    dat_in_filt(k) = fi(adc_input_int(j+k,z),1,12,0);
                    
                end 

                [y_out, y_out_int] = filter_transversal(dat_in_filt_double, www1, dat_in_filt, www1_int);

            end

            y_array(j+1,z-1) = y_out;
            y_array_int(j+1,z-1) = y_out_int;

            % if double(y_out_int)*2^-14 == -351
            %     eq = 1;
            %     figure(11);
            %     hold on
            %     subplot(2,1,1);
            %     plot(y_array(:,z-1));
            %     hold on
            %     subplot(2,1,2);
            %     plot(double(y_array_int(:,z-1)) *2^-14 );
            % 
            %     [qwe, qwe1] = filter_transversal(dat_in_filt_double, www1, dat_in_filt, www1_int);
            % end

        end


        % figure(12);
        % hold on
        % subplot(2,1,1);
        % plot([y_array(:,z-1)]);
        % hold on
        % subplot(2,1,2);
        % plot(double(y_array_int(:,z-1)) *2^-14 );


 end

 min_det_matrix = min(det_matlab);
 min_shift_det_matrix = min(det_x3_shift); 

 max_det_matrix = max(det_matlab);
 max_shift_det_matrix = max(det_x3_shift);

 minimum = min([min_det_matrix min_shift_det_matrix]);
 maximum = max([max_det_matrix max_shift_det_matrix]);

end