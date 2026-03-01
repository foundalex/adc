function [sig_adc, delta_tilda] = new_algorithm(adc_input, adc_input_sin, s_to_subadc_int, s_after_subadc, sim_options);


    % load ('2026              2             20             12              9          8.479.mat'); % f = 420MHz, ts = 0.05 
    % load ('2026              2             20             12             20         31.557.mat'); % f = 450 MHz, ts = 0.05

    N = 8192;
    % dt1 = 1/sim_options.Fs_sub_adc;
    hdc = 9;
    hdq1 = 25;
    hdq2 = 5;

    zf1 = zeros(sim_options.M,hdq1-1);
    zf2 = zeros(sim_options.M,hdq2-1);

    n = (0:1:sim_options.N-1);
    del_proc = ((sim_options.N-1)/2);
    ws1 = 0.998;
    ws2 = 0.001;
    w_blackman = 0.42 - 0.5 * cos(2*pi*n/(sim_options.N-1)) + 0.08 * cos(4*pi*n/(sim_options.N-1));

    lowpass_hilbert1 = (2./((n-del_proc)*pi)).*(sin(ws1.*((n-del_proc)*pi)./2)).^2;
    lowpass_hilbert1(1) = 0;
    lowpass_hilbert1(37) = 0;
    lowpass_hilbert2 = (2./((n-del_proc)*pi)).*(sin(ws2.*((n-del_proc)*pi)./2)).^2;
    lowpass_hilbert2(1) = 0;
    lowpass_hilbert2(37) = 0;

    weight_lowpass_hilbert1 = lowpass_hilbert1 .* w_blackman; 
    weight_lowpass_hilbert2 = lowpass_hilbert2 .* w_blackman; 
    bandpass_hilbert = weight_lowpass_hilbert1'  - weight_lowpass_hilbert2';

    %%
    bb = firpm(hdc,[0 0.9],[0 0.9*pi],'d');
    Hd = dfilt.dfasymfir(bb);

    bb1 = firpm(hdq1,[0.1 0.9],[1 1],'differentiator');
    Hd1(1) = dfilt.dfasymfir(bb1);
    % Hd1(2) = dfilt.dfasymfir(bb1);
    % Hd1(3) = dfilt.dfasymfir(bb1);
    % Hd1(4) = dfilt.dfasymfir(bb1);

    % Hd1(1,1).Numerator(12) = 0;
    % Hd1(1,2).Numerator(12) = 0;
    % Hd1(1,3).Numerator(12) = 0;
    % Hd1(1,4).Numerator(12) = 0;

    % figure(8);
    % subplot(1,1,1)
    % plot(Hd1(1,1).Numerator);

    Hd1s = dfilt.dfsymfir(bb1);

    % nn = 1:25;
    % figure(71); plot(nn, Hd1(1,1).Numerator);
    

    %%

    bb2 = firpm(hdq2,[0.1 0.9],[1 1],'differentiator');
    Hd2 = dfilt.dfasymfir(bb2);
    % blo = firls(18,[0 0.45 0.55 1],[1 1 0 0],[100 1], 'differentiator');

    [x, f] = freqz(Hd.Numerator, 1,1024, 'whole', sim_options.Fs_sub_adc);
    [x1, f1] = freqz(Hd1(1).Numerator, 1,1024, 'whole', sim_options.Fs_sub_adc);
    [xx1, ff1] = freqz(Hd1s.Numerator, 1,1024, 'whole', sim_options.Fs_sub_adc);
    
    [x2, f2] = freqz(Hd2.Numerator, 1,1024, 'whole', sim_options.Fs_sub_adc);

    % fvtool(Hd,Hd1,Hd2);

    %%
    % figure(3);
    % subplot(3,1,1)
    % plot(f, abs(x));
    % subplot(3,1,2)
    % plot(f1, abs(x1));
    % subplot(3,1,3)
    % plot(f2, abs(x2));
    % % 
    % figure(4);
    % subplot(4,1,1)
    % plot(f, angle(x));
    % subplot(4,1,2)
    % plot(f1, angle(x1));
    % title('ФЧХ фильтра дифф.фильтра 25-го порядка assym')
    % xlabel('Частота') 
    % ylabel('Амплитуда') 
    % subplot(4,1,3)
    % plot(ff1, angle(xx1));
    % title('ФЧХ фильтра дифф.фильтра 25-го порядка symm')
    % xlabel('Частота') 
    % ylabel('Амплитуда')
    % subplot(4,1,4)
    % plot(f2, angle(x2));
    %%
    % figure(5);
    % subplot(2,1,1)
    % plot(f1, angle(x1));
    % title('ФЧХ фильтра дифф.фильтра 25-го порядка')
    % xlabel('Частота') 
    % ylabel('Амплитуда') 
    % subplot(2,1,2)
    % plot(Hd1.Numerator);
    % title('Импульсная характеристика')
    % xlabel('Номер отсчета') 
    % ylabel('Амплитуда') 

    %%
    coeff(1,1) = 0;
    coeff(2,1) = 0;
    coeff(3,1) = 0;
    coeff(4,1) = 0;

    del_proc0 = (hdc-1)/2;
    del_proc1 = (hdq1-1)/2;
    del_proc2 = (hdq2-1)/2;

    U = [2 -1 0 -1; -1 2 -1 0; 0 -1 2 -1; -1 0 -1 2];
    % находим псевдоинверсную матрицу
	pseudo_u = pinv(U);

    % buffer = zeros(1,n0);

    y1_mult = zeros(N,sim_options.M);
    y2_mult = zeros(N,sim_options.M);
    window_adc = zeros(N,sim_options.M);
    window_signal = zeros(N,sim_options.M);
    %% 25 - ый порядок
    n1_impz = -(hdq1-1)/2:(hdq1-1)/2;
    n_hdq1 = 0:hdq1-1;

    w_blackman = 0.42 - 0.5 * cos(2*pi*n_hdq1/(hdq1-1)) + 0.08 * cos(4*pi*n_hdq1/(hdq1-1));
    for i = 1:hdq1
        hdq1_coeff(i) = ((-1)^n1_impz(i))/n1_impz(i); 
    end
    hdq1_coeff((hdq1+1)/2) = 0;

    % a1 == a2 
    % for i = 1:N_taps
    %     a2(i) = cos(n11(i)*pi)/n11(i); 
    % end
    % a2(1) = 0;
    % a2(13) = 0;

    % figure(71); 
    % plot([a1,a2]);

    diff_25 = hdq1_coeff .* w_blackman;

    Hdd(1) = dfilt.dfasymfir(diff_25);
    Hdd(2) = dfilt.dfasymfir(diff_25);
    Hdd(3) = dfilt.dfasymfir(diff_25);
    Hdd(4) = dfilt.dfasymfir(diff_25);

    % Hdd = dfilt.dffir(diff_25);
    % 
    [x3, f3] = freqz(diff_25, 1,1024, 'whole', sim_options.Fs_sub_adc);
    [x31, f31] = freqz(Hdd(1,1).Numerator, 1,1024, 'whole', sim_options.Fs_sub_adc);

    figure(70); 
    subplot(3,1,1)
    plot(hdq1_coeff)
    subplot(3,1,2)
    plot(f3, abs(x3), f3, abs(x31));
    subplot(3,1,3)
    plot(f3, angle(x3)); %, f3, angle(x31));

    %% 9 - ый порядок
    n_hdc_impz = -(hdc-1)/2:(hdc-1)/2;
    n_hdc = 0:hdc-1;

    w_blackman_hdc = 0.42 - 0.5 * cos(2*pi*n_hdc/(hdc-1)) + 0.08 * cos(4*pi*n_hdc/(hdc-1));
    for i = 1:hdc
        hdc_coeff(i) = ((-1)^n_hdc_impz(i))/n_hdc_impz(i); 
    end
    hdc_coeff((hdc+1)/2) = 0;

    diff_filter_hdc = hdc_coeff .* w_blackman_hdc;
   
    Hd_hdc = dfilt.dfasymfir(diff_filter_hdc);

    % 
    [x3_hdc, f3_hdc] = freqz(diff_filter_hdc, 1,1024, 'whole', sim_options.Fs_sub_adc);
    [x31_hdc, f31_hdc] = freqz(Hd_hdc.Numerator, 1,1024, 'whole', sim_options.Fs_sub_adc);

    figure(71); 
    subplot(3,1,1)
    plot(hdc_coeff)
    subplot(3,1,2)
    plot(f3_hdc, abs(x3_hdc), f3_hdc, abs(x31_hdc), f3_hdc, abs(x));
    subplot(3,1,3)
    plot(f3_hdc, angle(x3_hdc), f3_hdc, angle(x31_hdc), f3_hdc, angle(x));
    title('ФЧХ фильтра дифф.фильтра 9-го порядка')
    xlabel('Частота') 
    ylabel('Амплитуда')

    %% 5 - ый порядок
    n_hdq2_impz = -(hdq2-1)/2:(hdq2-1)/2;
    n_hdq2 = 0:hdq2-1;
    w_blackman_hdq2 = 0.42 - 0.5 * cos(2*pi*n_hdq2/(hdq2-1)) + 0.08 * cos(4*pi*n_hdq2/(hdq2-1));

    for i = 1:hdq2
        hdq2_coeff(i) = ((-1)^n_hdq2_impz(i))/n_hdq2_impz(i); 
    end
    hdq2_coeff((hdq2+1)/2) = 1;
    hdq2_coeff(4) = hdq2_coeff(2);
    hdq2_coeff(5) = hdq2_coeff(1);

    diff_filter_hdq2 = hdq2_coeff .* w_blackman_hdq2;
    Hd_hdq2 = dfilt.dfasymfir(diff_filter_hdq2);
    % 
    [x3_hdq2, f3_hdq2] = freqz(diff_filter_hdq2, 1,1024, 'whole', sim_options.Fs_sub_adc);
    [x31_hdq2, f31_hdq2] = freqz(Hd_hdq2.Numerator, 1,1024, 'whole', sim_options.Fs_sub_adc);

    figure(72); 
    subplot(3,1,1)
    plot(hdq2_coeff)
    subplot(3,1,2)
    plot(f3_hdq2, abs(x3_hdq2), f3_hdq2, abs(x31_hdq2), f3_hdq2, abs(x2));
    subplot(3,1,3)
    plot(f3_hdc, angle(x3_hdq2), f3_hdq2, angle(x31_hdq2), f3_hdq2, angle(x2));
    title('ФЧХ фильтра дифф.фильтра 5-го порядка')
    xlabel('Частота') 
    ylabel('Амплитуда')

    %%
    figure(80); 
    plot(f3_hdc, angle(x3_hdc), f3_hdc, angle(x3), f3_hdc, angle(x3_hdq2));


    %% VI samples intervention
    % in = 0;
    % row = 1;
    % start = 1;
    % 
    % while in < numel(s_after_subadc)-4
    %     for j = start:4
    %         in = in + 1;
    %         if in ~= 1 && mod(in-1,11) == 0
    %             j = 1;
    %             row = row + 1;
    %             sig(row,j) = s_after_subadc(in);
    %             start = 2;
    %         else
    %             sig(row,j) = s_after_subadc(in);
    %             start = 1;
    %         end
    %     end
    %     if mod(in-1,11) ~= 0
    %         row = row + 1;
    %     end
    % end
    % 
    % % удаление нулей из 4-го АЦП
    % odd_i = sig(1:3:end,4);
    % even_i = [sig(2:3:end,4)];
    % result = reshape([odd_i'; even_i'],1,[]);
    % 
    % adc_input_interv = zeros(numel(sig(:,1)), sim_options.M);
    % 
    % for i = 1:sim_options.M
    %     if i ~= sim_options.M
    %         adc_input_interv(:,i) = sig(:,i);
    %     else
    %         adc_input_interv(1:numel(result),i) = result;
    %     end
    % end

    %% H(d,c)
    % syms f(x)
    % f(x) = cos(x);
    % df = diff(f,x);

    for i = 1:sim_options.M
        y(:,i) = filter(Hd_hdc, double(adc_input(:,i)));
        % y(:,i) = filter(diff_filter_hdc, 1, double(adc_input(:,i)));
        y_gold(:,i) = filter([1 -1], 1, double(adc_input(:,i)));
        y_hilbert(:,i) = filter(bandpass_hilbert, 1, double(adc_input(:,i)));
        y_cut(:,i) = y(del_proc0+1:end,i);

        % figure(4)
        % subplot(3,1,1)
        % plot([double(adc_input(1:200,i)), y_cut(1:200,i)]);
        % subplot(3,1,2)
	    % snr(double(adc_input(:,i)), sim_options.Fs_sub_adc);
        % subplot(3,1,3)
	    % snr(y_cut(:,i), sim_options.Fs_sub_adc);

        % for n = 1:100
        %     r = r + 1;
        %     if r == 4
        %         arc_cos(n,i) = real(-acos(adc_input(n,i)));
        %         r = 0;
        %     else
        %         arc_cos(n,i) = real(acos(adc_input(n,i)));
        %     end
        %     % arc_cos(:,i) = real(acos(adc_input(n,i)));
        %     diff_func(n,i) = double(df(arc_cos(n,i)));
        % end
        % difference(:,i) = double(adc_input(1:100,i)) - input_s(1:100,i);
    end
    % 
    % figure(92); plot([adc_input(1:100,1), adc_input_sin(1:100,1)]);
    % figure(90); plot([adc_input(1:100,1), ...
    %     adc_input_sin(1:100,1), ...
    %     y_hilbert(37:136,1), ...
    %    ... diff_func(1:100,1), ...
    %     y_cut(1:100,1)]);
    % title('Сигналы ')
    % xlabel('Номер отсчета') 
    % ylabel('Амплитуда')
    % legend({'Входной сигнал cos(x)', 'Входной сигнал -sin(x) ', 'Выходной сигнал фильтра Гилберта', 'Выходной сигнал дифф.фильтра'}, 'Location','northwest');

    % [x3_gold, f3_gold] = freqz([1 -1], 1,1024, 'whole', sim_options.Fs_sub_adc);

    % figure(91);
    % subplot(2,1,1)
    % plot(f3_gold, abs(x3_gold), f3_hdc, abs(x3_hdc));
    % subplot(2,1,2)
    % plot(f3_gold, angle(x3_gold), f3_hdc, angle(x3_hdc));
    

    % syms f(x)
    % f(x) = cos(x);
    % df = diff(f,x);
    % eq = double(df(double(adc_input(10,1))*2^-11));
    % eq_int = eq * 2^11;

    % input_s = -asin(eq)*2^11;


    % load ('2026              2             23             23             49         25.385.mat'); % f = 420MHz, ts = 0.05 

    for j = 1:floor(length(double(y_cut(:,1)))/N) 

        start_index = N*(j-1)+1;
        end_index = start_index+N-1;

        % формируем окно в 8192 отсчета
        for i = 1:sim_options.M
            window_signal(:,i) = double(adc_input(start_index:end_index,i));
            window_adc(:,i) = double(y_cut(start_index:end_index,i));
        end

        % фильтруем
        for i = 1:sim_options.M

            y_mult(:,i) = window_adc(:,i) .* coeff(i,j);
            % y_mult(:,i) = window_adc(:,i) .* 1;

            y_tilda(:,i) = double(adc_input(start_index:end_index,i)) - y_mult(:,i); % разность входного сигнала и выхода первого фильтра y(n)-y'
            % y_tilda(:,i) = window_adc(:,i);

            % figure(7)
            % subplot(2,1,1)
            % plot(y_tilda(1:200));
            % subplot(2,1,2)
	        % snr(y_tilda(:,i), sim_options.Fs_sub_adc);

            %% 

            % [y1]  = filter(Hdd(i), y_tilda(:,i));

            % [y1, zf1(i,:)] = filter(diff_25, 1, y_tilda(:,i), zf1(i,:));
            % if j == 1
                % y1_cut(:,i) = [y1(del_proc1+1:end); zeros(del_proc1,1)];
            % else
                % y1_cut(:,i) = y1;
            % end

            % figure(5)
            % subplot(3,1,1)
            % plot([(y_tilda(1:200,i)), y1_cut(1:200,i)]);
            % subplot(3,1,2)
	        % snr((y_cut(:,i)), sim_options.Fs_sub_adc);
            % subplot(3,1,3)
	        % snr(y1_cut(:,i), sim_options.Fs_sub_adc);

            % y1_mult(:,i) = y1_cut(:,i) .* coeff(i,j);

            % y1_mult(:,i) = y1_cut(:,i) .* ceff1(i);
            y1_mult(:,i) =  y_tilda(:,i);


            % y2 = filter(Hd2, y_tilda(:,i));

            % [y2, zf2(i,:)] = filter(diff_filter_hdq2, 1, y_tilda(:,i), zf2(i,:));
            % if j == 1
            %     y2_cut(:,i) = [y2(del_proc2+1:end); zeros(del_proc2,1)];
            % else
            %     y2_cut(:,i) = y2;
            % end

            % figure(6)
            % subplot(3,1,1)
            % plot([(y_cut(1:200,i)), y2_cut(1:200)]);
            % subplot(3,1,2)
	        % snr((y_cut(:,i)), sim_options.Fs_sub_adc);
            % subplot(3,1,3)
	        % snr(y2_cut, sim_options.Fs_sub_adc);

            % y2_mult(:,i) = y2_cut(:,i) .* (coeff(i,j)^2)/2; 
            % y2_mult(:,i) = y2_cut(:,i) .* 1;

            x_tilda(start_index:end_index,i) = y1_mult(:,i); %double(adc_input(start_index:end_index,i)) - y1_mult(:,i); % - y2_mult(:,i);

            % figure(6)
            % subplot(2,1,1)
            % plot([adc_input(1:200,i), x_tilda(1:200,i)]);
            % subplot(2,1,2)
	        % snr(x_tilda(:,i), sim_options.Fs_sub_adc);

        end

        % производим корелляцию в окне и находим среднее значение
        cor1(:,j) = sum(x_tilda(start_index:end_index,1).*x_tilda(start_index:end_index,2))/N;
        cor2(:,j) = sum(x_tilda(start_index:end_index,2).*x_tilda(start_index:end_index,3))/N;
        cor3(:,j) = sum(x_tilda(start_index:end_index,3).*x_tilda(start_index:end_index,4))/N;
        cor4(:,j) = sum([0; x_tilda(start_index:end_index-1,4)] .* x_tilda(start_index:end_index,1))/N;

        % находим ошибку среднего значения корреляций
		er(1,j) = cor1(:,j) - cor4(:,j);
		er(2,j) = cor2(:,j) - cor1(:,j);
		er(3,j) = cor3(:,j) - cor2(:,j);
		er(4,j) = cor4(:,j) - cor3(:,j);


		%% выход АЦП0 дифф. фильтра

        cor_M(j) = sum(window_adc(:,1) .* x_tilda(start_index:end_index,2))/N;

		rate = nextpow2(cor_M(j));
		cor_M_floor = 2^rate;

		%%
		mult_u = pseudo_u * er(:,j);

		delta_tilda(:,j) = (mult_u/cor_M_floor);

        coeff(:,j+1) = coeff(:,j) + 0.2 * delta_tilda(:,j);


    end

    % save (sprintf(num2str(clock) + ".mat"));
    % load ('2026              2             23             23             51         10.257.mat'); % f = 420MHz, ts = 0.05 

    % выход АЦП
    % sig_adc = zeros(sim_options.M*length(x_tilda(:,1)),1);

    % for i = 1:numel(x_tilda(:,1))
    %     for j = start:4
    %         in = in + 1;
    %         if in ~= 1 && mod(in-1,11) == 0
    %             j = 1;
    %             row = row + 1;
    %             sig_adc(row) = x_tilda(i,j);
    %             start = 2;
    %         else
    %             row = row + 1;
    %             sig_adc(row) = x_tilda(i,j);
    %             start = 1;
    %         end
    %     end
    %     % if mod(in-1,11) ~= 0
    %     %     row = row + 1;
    %     % end
    % end

    % удаление нулей из 4-го АЦП

    % adc4 = x_tilda(1:363662,4);

    sig_adc = zeros(sim_options.M*length(x_tilda(:,1)),1);
    for i = 1:sim_options.M
        sig_adc(i:sim_options.M:end) = x_tilda(:,i);
    end

    % sig_adc = sig_adc(1:1452642);

    % result1 = reshape([odd_i'; even_i'],1,[]);

    figure(6)
    subplot(8,1,1)
    plot([s_to_subadc_int(end-500:end), s_after_subadc(end-500:end), sig_adc(end-500:end)]);
    title('Входной сигнал и выходной сигнал')
    xlabel('Номер отсчета') 
    ylabel('Значение отсчета') 
    legend({'Вход','Выход'},'Location','northeast');
    %% дельты
    subplot(8,1,2)
    plot(delta_tilda(1,:));
    title('Дельта первого суб-АЦП')
    xlabel('Номер итерации') 
    ylabel('Значение дельты') 
    subplot(8,1,3)
    plot(delta_tilda(2,:));
    title('Дельта второго суб-АЦП')
    subplot(8,1,4)
    plot(delta_tilda(3,:));
    title('Дельта третьего суб-АЦП')
    subplot(8,1,5)
    plot(delta_tilda(4,:));
    title('Дельта четвертого суб-АЦП')
    subplot(8,1,6)
    snr(s_to_subadc_int, sim_options.Fs_sub_adc*4);
    subplot(8,1,7)
    snr(s_after_subadc, sim_options.Fs_sub_adc*4);
    subplot(8,1,8)
	snr(sig_adc, sim_options.Fs_sub_adc*4);

    %%
    figure(7);
    subplot(4,1,1)
    plot(er(1,:));
    subplot(4,1,2)
    plot(er(2,:));
    subplot(4,1,3)
    plot(er(3,:));
    subplot(4,1,4)
    plot(er(4,:));

end