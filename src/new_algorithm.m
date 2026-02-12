function [] = new_algorithm(s_to_subadc_int, s_after_subadc, sim_options);

    N = 8192;
    n0 = 9;
    n1 = 25;
    n2 = 5;

    bb = firpm(n0,[0.1 0.9],[1 1],'differentiator');
    bb1 = firpm(n1,[0.1 0.9],[1 1],'differentiator');
    bb2 = firpm(n2,[0.1 0.9],[1 1],'differentiator');

    % blo = firls(18,[0 0.45 0.55 1],[1 1 0 0],[100 1], 'differentiator');

    [x, f] = freqz(bb, 1,1024, 'whole', sim_options.Fs_sub_adc);
    [x1, f1] = freqz(bb1, 1,1024, 'whole', sim_options.Fs_sub_adc);
    [x2, f2] = freqz(bb2, 1,1024, 'whole', sim_options.Fs_sub_adc);

    figure(3);
    subplot(3,1,1)
    plot(f, abs(x));
    subplot(3,1,2)
    plot(f1, abs(x1));
    subplot(3,1,3)
    plot(f2, abs(x2));

	% y = filter(bb, 1, double(s_to_subadc_int));
	% figure(4)
	% snr(y, sim_options.Fs_sub_adc*4);

    % yy = double(s_after_subadc) - y;
    % 
    % figure(4);
    % subplot(2,1,1)
    % plot([double(s_after_subadc(1:200)), y(1:200)]);
    % subplot(2,1,2)
    % plot([yy(1:200)]);
    % 
    % y1 = (filter(bb1, 1, yy));
    % 
    % y2 = (filter(bb2, 1, yy));
    % 
    % yyy = double(s_after_subadc) - y1 - y2;
    % 
	% figure(5);
	% subplot(3,1,1)
	% plot([y1(40:200)]);
	% subplot(3,1,2)
	% plot([yyy(40:200)]);
	% subplot(3,1,3)
	% plot([double(s_after_subadc(40:200))]);
    % 
	% figure(6)
	% snr(yyy, sim_options.Fs_sub_adc*4);

    % a1 = 1;

    % delta_tilda1 = ones(sim_options.M,1);
    delta_tilda1(1) = 0.8;
    delta_tilda1(2) = 0.8;
    delta_tilda1(3) = 0.8;
    delta_tilda1(4) = 0.8;

    del_proc0 = (n0-1)/2;
    del_proc1 = (n1-1)/2;
    del_proc2 = (n2-1)/2;

    for j = 1:floor(length(s_after_subadc)/N/sim_options.M) 

        adc = zeros(N,sim_options.M);
 
        start_index = sim_options.M*N*(j-1)+1;
        end_index = sim_options.M*N*(j-1)+N*sim_options.M;
        
        % формируем окно в 8192 отсчета
        for i = 1:sim_options.M
            adc(:,i) = s_after_subadc(start_index+i-1:sim_options.M:end_index);
        end

        % фильтруем
        for i = 1:sim_options.M
            y = filter(bb, 1, adc(:,i));
            y = y(del_proc0:end);
            y_mult(:,i) = y .* delta_tilda1(i);
            y_tilda(:,i) = adc(1:end-del_proc0+1,i) - y_mult(:,i); % разность входного сигнала и выхода первого фильтра y(n)-y'

            %%
            y1 = (filter(bb1, 1, y_tilda(:,i)));
            y1 = y1(del_proc1:end);

            y1_mult(:,i) = y1 .* delta_tilda1(i);

            y2 = (filter(bb2, 1, y_tilda(:,i)));
            y2 = y2(del_proc2:end);

            % y2_mult(:,i) = y2 .* delta_tilda1(i)/2;
            y2_mult(:,i) = y2 .* delta_tilda1(i); % исправить!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


            x_tilda(:,i) = adc(1:length(y1_mult(:,1)),i) - y1_mult(:,i) - y2_mult(1:length(y1_mult(:,1)),i);
        end


        % производим корелляцию в окне и находим среднее значение
        cor1 = sum(x_tilda(:,1).*x_tilda(:,2))/N;
		cor2 = sum(x_tilda(:,2).*x_tilda(:,3))/N;
		cor3 = sum(x_tilda(:,3).*x_tilda(:,4))/N;
        cor4 = sum([0; x_tilda(1:end-1,4)] .* x_tilda(:,1))/N;

        % находим ошибку среднего значения корреляций
		e(1,1) = cor1;
		e(2,1) = cor2 - cor1;
		e(3,1) = cor3 - cor2;
		e(4,1) = cor4 - cor3;

		%% 
        
        % выход АЦП0 дифф. фильтра
	
		cor_M = sum(y_tilda(1:length(x_tilda(:,1)),1) .* x_tilda(:,2))/N;
	
		rate = nextpow2(cor_M * -1);
		cor_M_floor = 2^rate;

		%%

		U = [2 -1 0 -1; -1 2 -1 0; 0 -1 2 -1; -1 0 -1 2];

        % находим псевдоинверсную матрицу
		pseudo_u = pinv(U);
	
		mult_u = pseudo_u * e;
		
		delta_tilda(:,j) = -mult_u/cor_M_floor;

		
    end


    sig_adc = zeros(sim_options.M*length(x_tilda(:,1)),1);
    for i = 1:sim_options.M
        sig_adc(i:sim_options.M:end) = x_tilda(:,i);
    end

    figure(6)
    subplot(3,1,1)
    plot([s_after_subadc(1:400), sig_adc(101:500)]);
    title('Входной сигнал и выходной сигнал')
    xlabel('Номер отсчета') 
    ylabel('Значение отсчета') 
    legend({'Вход','Выход'},'Location','northeast');
    subplot(3,1,2)
    snr(s_after_subadc, sim_options.Fs_sub_adc*4, 3);
    subplot(3,1,3)
	snr(sig_adc, sim_options.Fs_sub_adc*4, 3);



    % outputVmax = helperHarmonicDistortionAmplifier(sig_adc);
    % figure(60);
    % helperPlotPeriodogram(outputVmax, sim_options.Fs_sub_adc*4, "power","annotate");
end