function [sig_adc, delta_tilda] = new_algorithm(adc_input, s_to_subadc_int, s_after_subadc, sim_options);


    % load ('2026              2             20             12              9          8.479.mat'); % f = 420MHz, ts = 0.05 
    % load ('2026              2             20             12             20         31.557.mat'); % f = 450 MHz, ts = 0.05

    N = 8192;
    n0 = 9;
    n1 = 25;
    n2 = 5;

    zf1 = zeros(sim_options.M,n1);
    zf2 = zeros(sim_options.M,n2);
    States1 = zeros(sim_options.M,n1);

    bb = firpm(n0,[0.1 0.9],[1 1],'differentiator');
    Hd = dfilt.dfasymfir(bb);
    %%
    bb1 = firpm(n1-1,[0.1 0.95],[1 1],'differentiator');
    Hd1(1) = dfilt.dfasymfir(bb1);
    Hd1(2) = dfilt.dfasymfir(bb1);
    Hd1(3) = dfilt.dfasymfir(bb1);
    Hd1(4) = dfilt.dfasymfir(bb1);

    % Hd1(1,1).Numerator(12) = 0;
    % Hd1(1,2).Numerator(12) = 0;
    % Hd1(1,3).Numerator(12) = 0;
    % Hd1(1,4).Numerator(12) = 0;


    figure(8);
    subplot(1,1,1)
    plot(Hd1(1,1).Numerator);

    Hd1s = dfilt.dfsymfir(bb1);

    
    nn = 1:26;
    % figure(71); plot(nn, bb1, nn, Hd1.Numerator);

    %%

    bb2 = firpm(n2,[0.1 0.9],[1 1],'differentiator');
    Hd2 = dfilt.dfasymfir(bb2);
    % blo = firls(18,[0 0.45 0.55 1],[1 1 0 0],[100 1], 'differentiator');

    [x, f] = freqz(Hd.Numerator, 1,1024, 'whole', sim_options.Fs_sub_adc);
    [x1, f1] = freqz(Hd1(1).Numerator, 1,1024, 'whole', sim_options.Fs_sub_adc);
    [xx1, ff1] = freqz(Hd1s.Numerator, 1,1024, 'whole', sim_options.Fs_sub_adc);
    
    [x2, f2] = freqz(Hd2.Numerator, 1,1024, 'whole', sim_options.Fs_sub_adc);

    % fvtool(Hd,Hd1,Hd2);

    %%
    figure(3);
    subplot(3,1,1)
    plot(f, abs(x));
    subplot(3,1,2)
    plot(f1, abs(x1));
    subplot(3,1,3)
    plot(f2, abs(x2));
    % 
    figure(4);
    subplot(4,1,1)
    plot(f, angle(x));
    subplot(4,1,2)
    plot(f1, angle(x1));
    title('ФЧХ фильтра дифф.фильтра 25-го порядка assym')
    xlabel('Частота') 
    ylabel('Амплитуда') 
    subplot(4,1,3)
    plot(ff1, angle(xx1));
    title('ФЧХ фильтра дифф.фильтра 25-го порядка symm')
    xlabel('Частота') 
    ylabel('Амплитуда')
    subplot(4,1,4)
    plot(f2, angle(x2));

    %%
    coeff(1,1) = 0;
    coeff(2,1) = 0;
    coeff(3,1) = 0;
    coeff(4,1) = 0;

    del_proc0 = (n0-1)/2;
    del_proc1 = (n1-1)/2;
    del_proc2 = (n2-1)/2;

    U = [2 -1 0 -1; -1 2 -1 0; 0 -1 2 -1; -1 0 -1 2];
    % находим псевдоинверсную матрицу
	pseudo_u = pinv(U);

    % buffer = zeros(1,n0);

    y1_mult = zeros(N,sim_options.M);
    y2_mult = zeros(N,sim_options.M);
    window_adc = zeros(N,sim_options.M);

    %%
    N_taps = 25;
    n11 = 0:24;
    w_blackman = 0.42 - 0.5 * cos(2*pi*n11/(N_taps-1)) + 0.08 * cos(4*pi*n11/(N_taps-1));
    for i = 1:N_taps
        a1(i) = ((-1)^n11(i))/n11(i); 
    end
    a1(1) = 0;
    a1(13) = 0;

    %% a1 == a2
    % for i = 1:N_taps
    %     a2(i) = cos(n11(i)*pi)/n11(i); 
    % end
    % a2(1) = 0;
    % a2(13) = 0;

    % figure(71); 
    % plot([a1,a2]);
    %%

    diff_25 = a1 .* w_blackman;
    Hdd = dfilt.dfasymfir(diff_25);

    [x3, f3] = freqz(diff_25, 1,1024, 'whole', sim_options.Fs_sub_adc);
    [x31, f31] = freqz(Hdd.Numerator, 1,1024, 'whole', sim_options.Fs_sub_adc);

    figure(70); 
    subplot(3,1,1)
    plot(a1)
    subplot(3,1,2)
    plot(f3, abs(x3), f3, abs(x31));
    subplot(3,1,3)
    plot(f3, angle(x3), f3, angle(x31));

    %% H(d,c)
    % for i = 1:sim_options.M
    %     y(:,i) = filter(Hd, double(adc_input(:,i)));
    %     % y(:,i) = filter(bb, 1, double(adc_input(:,i)));
    % 
    %     y_cut(:,i) = y(del_proc0:end,i);
    % 
    %     % figure(4)
    %     % subplot(3,1,1)
    %     % plot([double(adc_input(1:200,i)), y_cut(1:200,i)]);
    %     % subplot(3,1,2)
	%     % snr(double(adc_input(:,i)), sim_options.Fs_sub_adc);
    %     % subplot(3,1,3)
	%     % snr(y_cut(:,i), sim_options.Fs_sub_adc);
    % end
    
    % ceff1 = [-0.0745, ...
    %             0.02, ...
    %             0.02, ...
    %             0.02];

     % ceff1 = [-0.274, ...
     %            0.09, ...
     %            0.09, ...
     %            0.09];
       % 
       % ceff1 = [-0.02, ...
       %          0.00, ...
       %          0.00, ...
       %          0.00];

    % %% VI samples intervention
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
    %%

    % load ('2026              2             23             23             49         25.385.mat'); % f = 420MHz, ts = 0.05 

    for j = 1:floor(length(adc_input(:,1))/N) 

        start_index = N*(j-1)+1;
        end_index = start_index+N-1;

        % формируем окно в 8192 отсчета
        for i = 1:sim_options.M
            window_adc(:,i) = double(adc_input(start_index:end_index,i));
        end

        % фильтруем
        for i = 1:sim_options.M

            % y_mult(:,i) = window_adc(:,i) .* coeff(i,j);
            % y_mult(:,i) = window_adc(:,i) .* 1;

            % y_tilda(:,i) = double(adc_input(start_index:end_index,i)) - y_mult(:,i); % разность входного сигнала и выхода первого фильтра y(n)-y'
            y_tilda(:,i) = window_adc(:,i);

            % figure(7)
            % subplot(2,1,1)
            % plot(y_tilda(1:200));
            % subplot(2,1,2)
	        % snr(y_tilda(:,i), sim_options.Fs_sub_adc);

            %% 

            % y1 = filter(Hd1, y_tilda(:,i));

   
            [y1]  = filter(Hd1(i), y_tilda(:,i));
            % Hd1.States = States1(:,i);
            % States1(:,i) = Hd1.States;
            % zerophase(y1, 1);

            % [y1, zf1(i,:)] = filter(bb1, 1, y_tilda(:,i), zf1(i,:));
            % if j == 1
                % y1_cut(:,i) = [y1(del_proc1:end); zeros(del_proc1-1,1)];
            % else
                y1_cut(:,i) = y1;
            % end

            % figure(5)
            % subplot(3,1,1)
            % plot([(y_cut(1:200,i)), y1_cut(1:200,i)]);
            % subplot(3,1,2)
	        % snr((y_cut(:,i)), sim_options.Fs_sub_adc);
            % subplot(3,1,3)
	        % snr(y1_cut(:,i), sim_options.Fs_sub_adc);

            y1_mult(:,i) = y1_cut(:,i) .* coeff(i,j);

            % y1_mult(:,i) = y1_cut(:,i) .* ceff1(i);
            % y1_mult(:,i) = y1_cut(:,i) .* 0;


            % y2 = filter(Hd2, y_tilda(:,i));

            % [y2, zf2(i,:)] = filter(bb2, 1, y_tilda(:,i), zf2(i,:));
            % if j == 1
                % y2_cut(:,i) = [y2(del_proc2:end); zeros(del_proc2-1,1)];
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

            x_tilda(start_index:end_index,i) = double(adc_input(start_index:end_index,i)) - y1_mult(:,i); % - y2_mult(:,i);

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

        cor_M = sum(y1_cut(:,1) .* x_tilda(start_index:end_index,2))/N;

		rate = nextpow2(cor_M);
		cor_M_floor = 2^rate;

		%%
		mult_u = pseudo_u * er(:,j);

		delta_tilda(:,j) = -(mult_u/cor_M_floor);

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
    %%
    subplot(8,1,6)
    snr(s_to_subadc_int, sim_options.Fs_sub_adc*4);
    subplot(8,1,7)
    snr(s_after_subadc, sim_options.Fs_sub_adc*4);
    subplot(8,1,8)
	snr(sig_adc, sim_options.Fs_sub_adc*4);


end