function yri_cut = fractional_delays(input_signal, input_signal_int, M, N_taps, Z)

    n = (0:1:N_taps-1);
    Nbp = floor(Z/2); % стр 7. (24)
    nn = 1:length(input_signal(:,1));
    del_proc = ((N_taps-1)/2);
    delay_adc = (1/M:1/M:1); % (стр.6,(16)), создаем массив на различные значения задержек

    w_blackman = 0.42 - 0.5 * cos(2*pi*n/(N_taps-1)) + 0.08 * cos(4*pi*n/(N_taps-1)); % Blackman window


    input_signal_fix = fi(input_signal, 1, 12, 10);
    % bb = double(input_signal_fix);
    % bb1 = double(input_signal_int) * 2^-11;
    % 
    % figure(3);
    % plot([input_signal(1:100), bb(1:100), bb1(1:100)])




    for i = 1:M-1
        D = del_proc - delay_adc(i); % delay (N-1)/2 - d = causal filter

        w_blackman_fractional(:,i) = 0.42 - 0.5 * cos(2*pi*(n+delay_adc(i))/(N_taps-1)) + 0.08 * cos(4*pi*(n+delay_adc(i))/(N_taps-1)); % shift Blackman window

        %% int
        % unsigned 8 bit
        w_blackman_fractional_int(:,i) = int16(round(w_blackman_fractional(:,i) * 2^8)); % fix(0,8,8)
        aa = double(w_blackman_fractional_int(:,i))*2^-8;

        % unsigned 8 bit
        w_blackman_fractional_fix(:,i) = fi(w_blackman_fractional(:,i),0,8,8); % fix(0,8,8)
        aa1 = double(w_blackman_fractional_fix(:,i));

        % figure(2);
        % plot([w_blackman_fractional(:,i), aa, aa1]);
        %%
        hri_m(:,i) = sinc(n-D); % shift impulse response on D = Dint - d for fractional delay filter

        %%
        % fix(1,9,8)
        hri_fix(:,i) = fi(hri_m(:,i),1,9,8);
        aa2 = double(hri_fix(:,i));
        % signed 9 bit
        hri_int(:,i) = int16(round(hri_m(:,i) * 2^9));
        aa3 = double(hri_int(:,i))*2^-9;

        % figure(4);
        % plot([hri_m(:,i), aa3, aa2]);
        %%
         
        hri_m(:,i) = hri_m(:,i) .* w_blackman_fractional(:,i); 
        % hri_m(:,i) = hri_m(:,i) ./ sum(hri_m(:,i));

        %%
        % fi(1,17,0)
        hri_mult_fix(:,i) = hri_fix(:,i) .* w_blackman_fractional_fix(:,i);
        hri_mult_int(:,i) = int32(hri_int(:,i)) .* int32(w_blackman_fractional_int(:,i));

        aa4 = double(hri_mult_int(:,i))*2^-17;
        aa5 = double(hri_mult_fix(:,i));

        % figure(4);
        % plot([hri_m(:,i), aa4 aa5,]);
        %%
        [yy(:,i), ff] = freqz(hri_m(:,i),1,1024, 'whole', 1000000000);

        yri(:,i) = filter(hri_m(:,i), 1, input_signal(:,1)); % (стр.6 (15))

        y_out_array_fix(:,i) = filter(hri_mult_fix(:,i), 1, input_signal_fix(:,1)); % (стр.6 (15))
        y_out_array_int(:,i) = int32(filter(hri_mult_int(:,i), 1, input_signal_int(:,1))); % (стр.6 (15))

        y_out_array_fix(:,i) = fi(y_out_array_fix(:,i),1,16,11);

        % 16+12
        y_out_array_int16(:,i) = (int16(round(y_out_array_int(:,i)/65536))).'; % fi(1,16,12)
        %%
        % figure(4);
        % plot([yri(1:200,i), double(y_out_array_fix(1:200,i)), double(y_out_array_int16(1:200,i)) * 2^-12]);

        figure(3);
        subplot(3,1,1);
        sfdr(yri(:,i), 1000000000);
        subplot(3,1,2);
        sfdr(double(y_out_array_fix(:,i)).', 1000000000);
        subplot(3,1,3);
        sfdr((double(y_out_array_int16(:,i))) * 2^-12, 1000000000);

        %% Thiran IIR All-pass
        % sys = thiran(D,1);
        % freqz(cell2mat(sys.Numerator),cell2mat(sys.Denominator),1024, 'whole', 1000000000);
        % yri(:,i) = filter(cell2mat(sys.Numerator), cell2mat(sys.Denominator), input_signal(:,1));
        % grpdelay(cell2mat(sys.Numerator), cell2mat(sys.Denominator),256,'whole', 1000000000);
      
        % h250 = h250 .* w_blackman_fractional.'; 
        % % h250 = normalize(h250,"norm",2);
        % h250 = h250/2.5;
     
        % % [aa, ff1] = freqz(hri_m(:,i),1,1024, 'whole', 1000000000);
        % % [aa1, ff2] = freqz(h250,1,1024, 'whole', 1000000000);
        % plot(ff1, abs(aa), ff1, abs(aa1))
        % % grpdelay(hri_m(:,i),1,256,'whole', 1000000000);
        % % grpdelay(h250,1,256,'whole', 1000000000);
        % 
        % aba(:,i) = filter(h250, 1, input_signal(:,1)); % (стр.6 (15))
        % 
        % yri(:,i) = real(aba(:,i));
        % plot([input_signal(1:500,2), (yri(1:500))]);

        %% Algorithm for working in different zones of Nyquist

        hh = (2./((n-del_proc)*pi)).*(sin(((n-del_proc)*pi)./2)).^2;
        hh(1) = 0;
        hh(37) = 0;

        hh_m = (hh .* w_blackman).';

        hh_m_fix = fi(hh_m, 1,8,7);

        % signed 8 bit
        hh_m_int = int16(round(hh_m * 2^7)); 
        % hh_m_int1 = double(hh_m_int) * 2^-8;

        % figure(4);
        % plot([hh_m,  double(hh_m_fix), hh_m_int]);


        ymi = filter(hh_m.', 1, yri(:,i));
        ymi_int = int32((filter(hh_m_int, 1, y_out_array_int16(:,i)))); % (стр.6 (15))

        ymi_fix = filter(hh_m_fix, 1, y_out_array_fix(:,i)); % (стр.6 (15))
        ymi_fix = fi(ymi_fix,1,14,12);

        % fi(1,16,12) + fi(1,8,7)
        ymi_d = double(ymi_int) * 2^-19;

        % 128 - 7 bit. fi(1,21,19) = fi(1,14,12)
        ymi_int16 = int16((floor(ymi_int/128))); 

        % figure(4);
        % plot([ymi(1:200), ymi_fix(1:200), ymi_d(1:200), double(ymi_int16(1:200) * 2^-12)]);

        figure(3);
        subplot(4,1,1);
        sfdr(ymi, 1000000000);
        subplot(4,1,2);
        sfdr(double(ymi_fix), 1000000000);
        subplot(4,1,3);
        sfdr(ymi_d, 1000000000);
        subplot(4,1,4);
        sfdr(double(ymi_int16) * 2^-12, 1000000000);

        %%
        yhil_imag(:,i) = [ymi(del_proc+1:end); zeros(del_proc,1)];
        yhil_imag_int(:,i) = [ymi_int16(del_proc+1:end); zeros(del_proc,1)];

        nn1 = nn + delay_adc(i);
        cc = cos(2*pi*nn1*Nbp);
        ss = sin(2*pi*nn1*Nbp);

        %%
        cc_fix = fi(cc,1,8,6);
        ss_fix = fi(ss,1,8,6);

        % figure(4);
        % plot([cc(1:200), double(cc_fix(1:200)), ss(1:200), double(ss_fix(1:200))]);
        %%
        cc_int = int8(round(cc * 2^6));
        ss_int = int8(round(ss * 2^6));
        %%
        yri_cos = cc.' .* yri(:,i); 
        yri_sin = ss.' .* yhil_imag(:,i); 
        %% 
        yri_cos_int = int32(cc_int.') .* int32(y_out_array_int16(:,i)); % fi(1,8,6) + fi(1,16,12) = fi(1,24,18);
        yri_cos_int16 = int16(round(yri_cos_int/256)); % fi(1,24,18) = fi(1,16,10)
        % figure(4);
        % plot([yri_cos(1:200), double(yri_cos_int16(1:200))*2^-10]);

        yri_sin_int = int32(ss_int.') .* int32(yhil_imag_int(:,i)); % fi(1,8,6) + fi(1,14,12) = fi(1,22,18);
        yri_sin_int16 = int16(round(yri_sin_int/256)); % fi(1,22,18) = fi(1,14,10)

        % figure(4);
        % plot([yri_sin(1:200), double(yri_sin_int16(1:200)*2^-10)]);
        % 
        % figure(3);
        % subplot(4,1,1);
        % sfdr(yri_cos, 1000000000);
        % subplot(4,1,2);
        % sfdr(yri_sin, 1000000000);
        % subplot(4,1,3);
        % sfdr(double(yri_cos_int16)*2^-10, 1000000000);
        % subplot(4,1,4);
        % sfdr(double(yri_sin_int16)*2^-10, 1000000000);

        %%
        if (mod(Z,2) == 0)
            yric(:,i) = yri_cos + yri_sin;
            yric_int(:,i) = yri_cos_int16 + yri_sin_int16; % fi(1,16,10) + fi(1,14,10) = fi(1,17,11)
        else
            yric(:,i) = yri_cos - yri_sin;
            yric_int(:,i) = yri_cos_int16 - yri_sin_int16; % fi(1,16,10) - fi(1,14,10) = fi(1,16,10)

            yric_int16 = double(yric_int)*2^-10;
            % figure(4);
            % plot([yric(1:200), yric_int16(1:200)]);

            figure(3);
            subplot(2,1,1);
            sfdr(yric, 1000000000);
            subplot(2,1,2);
            sfdr(yric_int16, 1000000000);
        end

        yri_cut(:,i+1) = yric(del_proc+1:end,i);
        yri_cut_int(:,i+1) = yric_int(del_proc+1:end,i);
    
    end
    yri_cut(:,1) = input_signal(1:end-del_proc,1);
    yri_cut_int(:,1) = input_signal_int(1:end-del_proc,1);

    %% example
    % hri_m1 = w_blackman .* sinc(n-del_proc); % shift impulse response on Dint
    % plot(n.',w_blackman.', n.',w_blackman_fractional(:,4).', n.',hri_m(:,4), n.',hri_m1)
    % legend({'Окно Блэкмена', 'Окно Блэкмена, сдвинутое на -0.5', 'Импульсная х-ка фильтра -0.5', 'Импульсная х-ка фильтра'},'Location','northeast')
    % xlabel('Номер отсчета') 
    % ylabel('Амплитуда') 

    %%
    % figure(8)
    % plot(ff, abs(yy(:,1)), ff, abs(yy(:,2)), ff, abs(yy(:,3)), ff, abs(yy(:,4)), ff, abs(yy(:,5)), ff, abs(yy(:,6)), ff, abs(yy(:,7)));
    % legend({'0.125','0.25', '0.375', '0.5', '0.625', '0.75', '0.875'},'Location','northeast')
    % title('АЧХ фильтров дробной задержки')
    % xlabel('Частота, Гц') 
    % ylabel('Коэффициент передачи') 
    % x81 = xline(5*10^8, '--', 'Fs/2');
    % x81.LabelHorizontalAlignment = 'center'
    % x81.LabelVerticalAlignment = 'middle';
    % x82 = xline(0.99*10^9, '--', 'Fs');
    % x82.LabelHorizontalAlignment = 'center'
    % x82.LabelVerticalAlignment = 'middle';

    % figure(9);
    % subplot(2,1,1)
    % plot([[zeros(del_proc,1); input_signal(1:274,2)], yric(1:del_proc+274,1)]);
    % % subplot(3,1,2)
    % % plot([input_signal(1:274,1), yri_cut(1:274,1)]);
    % subplot(2,1,2)
    % plot([input_signal(1:274,2), yri_cut(1:274,2)]);
    % subplot(9,1,4)
    % plot([input_signal(1:274,3), yri_cut(1:274,3)]);
    % subplot(9,1,5)
    % plot([input_signal(1:274,4), yri_cut(1:274,4)]);
    % subplot(9,1,6)
    % plot([input_signal(1:274,5), yri_cut(1:274,5)]);
    % subplot(9,1,7)
    % plot([input_signal(1:274,6), yri_cut(1:274,6)]);
    % subplot(9,1,8)
    % plot([input_signal(1:274,7), yri_cut(1:274,7)]);
    % subplot(9,1,9)
    % plot([input_signal(1:274,8), yri_cut(1:274,8)]);

    % spectrumScope = spectrumAnalyzer(SampleRate=1000000000, ...            
    %             AveragingMethod='exponential', ForgettingFactor=0, ...
    %             YLimits=[-30 10],ShowLegend=true, Method='Welch');
    % spectrumScope.WindowLength = 2048;
    % spectrumScope.FrequencyResolutionMethod = "window-length";
    % spectrumScope.PlotAsTwoSidedSpectrum=true;
    % spectrumScope.DistortionMeasurements.Enabled = true;
    % 
    % spectrumScope([input_signal(100:4096+100,2), yric(100:4096+100,1)]);

    %% test signal after fractional delay filters

    yri_cut1(:,1) = double(yri_cut_int(:,1))*2^-11; 
    for j = 2:M
        yri_cut1(:,j) = double(yri_cut_int(:,j)) * 2^-10;
    end

    % figure(9);
    % subplot(8,1,1)
    % plot([yri_cut(1:200,1), yri_cut1(1:200,1)]);
    % subplot(8,1,2)
    % plot([yri_cut(1:200,2), yri_cut1(1:200,2)]);
    % subplot(8,1,3)
    % plot([yri_cut(1:200,3), yri_cut1(1:200,3)]);
    % subplot(8,1,4)
    % plot([yri_cut(1:200,4), yri_cut1(1:200,4)]);
    % subplot(8,1,5)
    % plot([yri_cut(1:200,5), yri_cut1(1:200,5)]);
    % subplot(8,1,6)
    % plot([yri_cut(1:200,6), yri_cut1(1:200,6)]);
    % subplot(8,1,7)
    % plot([yri_cut(1:200,7), yri_cut1(1:200,7)]);
    % subplot(8,1,8)
    % plot([yri_cut(1:200,8), yri_cut1(1:200,8)]);

	sig_adc = zeros(M*length(yri_cut(:,1)),1);
    sig_adc_int = zeros(M*length(yri_cut_int(:,1)),1);
    
	for i = 1:M
		sig_adc(i:M:end) = yri_cut(:,i);
        sig_adc_int(i:M:end) = yri_cut1(:,i);
    end

    figure(4);
    plot([sig_adc(1:500), sig_adc_int(1:500)]);

    figure(3);
    subplot(2,1,1);
    sfdr(sig_adc, 1000000000);
    subplot(2,1,2);
    sfdr(sig_adc_int, 1000000000);



end