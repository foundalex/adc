function [y_array, y_array_int] = least_mean_square(adc_input, adc_input_int, yri_cut, yri_cut_int, M, N, width)

kk = 0;
tt = 0;

length_word_div = 70;

 for z = 2:M
    %% блок для расчета первых N коэффициентов фильтра
    %%
    
    % создаем матрицу входного сигнала
    for i = 1:N
        x3(i,:) = adc_input(i:N+i-1,z).'; % (стр.6, (20))
        x3_int(i,:) = adc_input_int(i:N+i-1,z).'; % fi(1,12,11)
    end

    % рассчитываем первые N коэффициентов адаптивного фильтра
    % сравнивая с задержанным сигналом ADC0 (yri_cut)
    % w1 = (x3' * x3) \ x3' * yri_cut(1:N,z); % (стр.6, (19))
    % w1 = lsqminnorm(x3, yri_cut(1:N,z));
   
    tt = tt + 1;

    %% initial determinant
    det_matlab(tt) = det(x3);
    [det_x3, det_x3_int] = determinate(x3, x3_int);

    % if (det_x3_int == 0)
    %     det_x3_int = fi(1,1,n,0);
    % end

    if (det_matlab(tt) == 0)
        det_matlab(tt) = 1;
    end

        %%
        for i = 1:N
            kk = kk + 1;

            x3_shift = x3;
            x3_shift(1:N,i) = yri_cut(1:N,z);

            % int
            x3_shift_int = x3_int;
            x3_shift_int(1:N,i) = yri_cut_int(1:N,z); 


            %% determinant
            det_x3_shift(kk) = det(x3_shift);
            [det_out_shift, det_out_shift_int] = determinate(x3_shift, x3_shift_int);

            if (det_x3_shift(kk) == 0)
                det_x3_shift(kk) = 1;
            end

            %% divide determinant
            www1(:,i) = det_x3_shift(kk) / det_matlab(tt); % double
            www1_int(:,i) = int_division(det_out_shift_int, det_x3_int, length_word_div, width); % integer

            www1_int_double(i,:) = double(www1_int(:,i))*2^-width;
        end

        %% filter
        % умножаем входные слова на рассчитанные коэффициенты
        % y_out = w1(1)*x(j+1) + w1(2)*x(j+2) + w1(3)*x(j+3) +  w1(4)*x(j+4); % Behrouz Farhang-Boroujeny, Adaptive Filters Theory and Applications  (стр. 414)

        for k = 1:N
            dat_in_filt_double(k) = adc_input(k,z);
            dat_in_filt(k) = fi(adc_input_int(k,z),1,12,0);
        end

        [y_out, y_out_int] = filter_transversal(dat_in_filt_double, www1, dat_in_filt, www1_int, width);

        %% array out
        y_array(1,z-1) = y_out;
        y_array_int(1,z-1) = y_out_int;




        %% part 2
        %%

        for j = 1:length(yri_cut(:,1))-2*N
  
            % shift to left matrix input signal. Refresh matrix input signal for every new word
            for i = 1:N
                x3(i,:) = [x3(i,2:N), 0];
                x3(i,N) = adc_input(j+N-1+i,z); % (стр.6, (20))

                x3_int(i,:) = [x3_int(i,2:N), 0]; % integer
                x3_int(i,N) = adc_input_int(j+N-1+i,z); 
            end
            %% determinant
           
            tt = tt + 1;
            det_matlab(tt) = det(x3);
            [det_x3(tt), det_x3_int(tt)] = determinate(double(x3_int)*2^-11, x3_int); % int

            if (det_matlab(tt) == 0)
                det_matlab(tt) = 1;
            end

            % if (j==17 & i ==5 & tt ==18)
            %     eq = 1;
            % end

            for i = 1:N
                kk = kk + 1;

                x3_shift = x3; % double
                x3_shift(1:N,i) = yri_cut(j+1:N+j,z); 

                x3_shift_int = x3_int; % int
                x3_shift_int(1:N,i) = yri_cut_int(j+1:N+j,z); 

                %% determinant

                det_x3_shift(kk) = det(x3_shift);
                [det_out_shift(kk), det_out_shift_int(kk)] = determinate(double(x3_shift_int)*2^-11, x3_shift_int);

                if (det_x3_shift(kk) == 0)
                    det_x3_shift(kk) = 1;
                end

                %% divide determinant
                % if (det_out_shift_int(kk) == 0)
                %     det_out_shift_int(kk) = fi(1,1,n,0);
                % end

                % if (j == 346 & i == 5 & z == 3)
                %     w = 1;
                % end  
                www1(:,i) = det_x3_shift(kk) ./ det_matlab(tt); % double

                www1_int(:,i) = int_division(det_out_shift_int(kk), det_x3_int(tt), length_word_div, width); % int
                www1_int_double(i,:) = double(www1_int(:,i))*2^-width;

                %% filter

                % y_outd = 0;

                % filter input signal. Mult input words on coeff
                for k = 1:N
                    % y_outd = y_outd + www1(k) * adc_input(j+k,z); % (стр 5, (13))

                    dat_in_filt_double(k) = adc_input(j+k,z);
                    dat_in_filt(k) = fi(adc_input_int(j+k,z),1,12,0);
                    
                end 

                [y_out, y_out_int] = filter_transversal(dat_in_filt_double, www1, dat_in_filt, www1_int, width);

            end

            y_array(j+1,z-1) = y_out;
            y_array_int(j+1,z-1) = y_out_int;

            % % % % if abs(y_out) > 100
            % if  abs(double(y_out_int)*2^0) > 2500
            %     eq = 1;
            %     figure(11);
            %     subplot(2,1,1);
            %     plot(y_array(:,z-1));
            %     subplot(2,1,2);
            %     plot(double(y_array_int(:,z-1)) *2^0);
            % 
            %     [qwe, qwe1] = filter_transversal(dat_in_filt_double, www1, dat_in_filt, www1_int);
            % end

        end
     
        % figure(13);
        % 
        % subplot(4,1,1);
        % plot(double(yri_cut_int(:,1)));
        % title('Выход референсного канала АЦП')
        % xlabel('Номер отсчета') 
        % ylabel('Амплитуда сигнала') 
        % 
        % subplot(4,1,2);
        % plot([y_array(:,z-1)]);
        % title('Выход адаптивного фильтра double для одного канала АЦП')
        % xlabel('Номер отсчета') 
        % ylabel('Амплитуда сигнала') 
        % 
        % subplot(4,1,3);
        % plot(double(adc_input_int(:,z)));
        % title('Вход для одного канала АЦП')
        % xlabel('Номер отсчета') 
        % ylabel('Амплитуда сигнала') 
        % 
        % subplot(4,1,4);
        % plot(double(y_array_int(:,z-1)) *2^0);
        % title('Выход адаптивного фильтра integer для одного канала АЦП')
        % xlabel('Номер отсчета') 
        % ylabel('Амплитуда сигнала') 

 end

        % figure(13);
        % subplot(4,1,1);
        % plot(double(yri_cut_int(:,1)));
        % subplot(4,1,2);
        % plot(double(y_array_int(:,1)));
        % subplot(4,1,3);
        % plot(double(y_array_int(:,2)));
        % subplot(4,1,4);
        % plot(double(y_array_int(:,3)));

 min_det_matrix = min(det_matlab);
 min_shift_det_matrix = min(det_x3_shift); 

 max_det_matrix = max(det_matlab);
 max_shift_det_matrix = max(det_x3_shift);

 minimum = min([min_det_matrix min_shift_det_matrix]);
 maximum = max([max_det_matrix max_shift_det_matrix]);

end