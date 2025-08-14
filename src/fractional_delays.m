function yri_cut = fractional_delays(input_signal, M, N_taps, Z)

    n = (0:1:N_taps-1);
    Nbp = floor(Z/2); % стр 7. (24)
    nn = 1:length(input_signal(:,1));
    del_proc = ((N_taps-1)/2);
    delay_adc = (1/M:1/M:1); % (стр.6,(16)), создаем массив на различные значения задержек
    w_blackman = 0.42 - 0.5 * cos(2*pi*n/(N_taps-1)) + 0.08 * cos(4*pi*n/(N_taps-1)); % Blackman window

    for i = 1:M-1
        D = del_proc - delay_adc(i); % delay (N-1)/2 - d = causal filter

        w_blackman_fractional(:,i) = 0.42 - 0.5 * cos(2*pi*(n+delay_adc(i))/(N_taps-1)) + 0.08 * cos(4*pi*(n+delay_adc(i))/(N_taps-1)); % shift Blackman window
        
        hri_m(:,i) = sinc(n-D); % shift impulse response on D = Dint - d for fractional delay filter
        hri_m(:,i) = hri_m(:,i) .* w_blackman_fractional(:,i); 
        [yy(:,i), ff] = freqz(hri_m(:,i),1,1024, 'whole', 1000000000);

        % hri_m(:,i) = hri_m(:,i) ./ sum(hri_m(:,i)); 
        yri(:,i) = filter(hri_m(:,i), 1, input_signal(:,1)); % (стр.6 (15))

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
        % yhil = hilbert(yri(:,i)); 
        % yhil_imag(:,i) = imag(yhil);

        hh = (2./((n-del_proc)*pi)).*(sin(((n-del_proc)*pi)./2)).^2;
        hh(1) = 0;
        hh(37) = 0;
        hh_m = hh .* w_blackman;
        ymi = filter(hh_m.', 1, yri(:,i));

        yhil_imag(:,i) = [ymi(del_proc+1:end); zeros(del_proc,1)];

        nn1 = nn + delay_adc(i);

        cc = cos(2*pi*nn1*Nbp);
        ss = sin(2*pi*nn1*Nbp);

        yri_cos = cc.' .* yri(:,i); 
        yri_sin = ss.' .* yhil_imag(:,i); 

        if (mod(Z,2) == 0)
            yric(:,i) = yri_cos + yri_sin;
        else
            yric(:,i) = yri_cos - yri_sin;
        end

        yri_cut(:,i+1) = yric(del_proc+1:end,i);
    
    end
    yri_cut(:,1) = input_signal(1:end-del_proc,1);
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
end