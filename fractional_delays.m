function [yri_cut, te] = fractional_delays(input_signal, M, N_taps)

    % создаем массив [1:3] из задержанного сигнала ADC0 на различные значения [0.75 0.5 0.25]
    % for i = 1:M-1
    %     delay_adc(i) = i/M; % (стр.6,(16))
    % end
    % delay_adc = fliplr(delay_adc);

    n = (1:N_taps);
    ww = blackman(N_taps);
    % % wvtool(blackman(N_taps));
    % 
    % % Nbp = 0; % if 1 zone of Nyquist
    % Nbp = 1; % if 2 zone of Nyquist
    % nn = 1:5000;

    for i = 1:M-1
        % Method 1
        delay_adc(i) = i/M; % (стр.6,(16))

        hri = sinc(n - (N_taps - 1) / 2 - delay_adc(i));
        te = hri.' .* ww;
        te = te ./ sum(te); 
        % freqz(te,1);
        yri = filter(te, 1, input_signal); % (стр.6 (15))

        % Method 2
        % [hri, i0, bw] = designFracDelayFIR(delay_adc(i), N_taps); 
        % [H1,w] = freqz(hri,1);
        % plot(w/pi,mag2db(abs([H1])))
        % total_delay_fir(i) = i0(i) + delay_adc(i);

        % plot([yri(1:100,i), adc_input(1:100,1)]);

        %% Algorithm for working in different zones of Nyquist
        % plot([adc_input(1:70,1), yri(1:70,i)]);
        % nn1 = nn + delay_adc(i);
        % cc(:,i) = cos(2*pi*Nbp*nn1);
        % ss(:,i) = sin(2*pi*Nbp*nn1);
        % y_hil(:,i) = hilbert(yri(:,i));
        % plot([real(y_hil(1:100)), imag(y_hil(1:100))]);
        % yri_cos = cc(:,i) .* yri(:,i); 
        % % plot(yri_cos(1:100))
        % yri_sin = ss(:,i) .* y_hil(:,i); 
        % sum_yri(:,i) = yri_cos + yri_sin;
        % % plot([adc_input(1:150,2), yri(1:150,i)]);
        % zz = sum_yri((N-1)/2:end,i);
        % % adc_input_id(1:end-(N-1)/2+1,i) = zz;
        % adc_input_id(:,i) = real(zz);

        % yri_cut = yri;
        yri_cut = yri((N_taps-1)/2:end);






        % yri_cut(:,i) = yri((N_taps-1)/2:end); % удаляем переходные процессы (N-1)/2
        % plot([adc_input(1:100,2), zz(1:100)]);

    end
end