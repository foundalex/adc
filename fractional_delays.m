function [yri_cut, te] = fractional_delays(input_signal, M, N_taps, Z)

    n = (1:N_taps);
    ww = blackman(N_taps);
    % % wvtool(blackman(N_taps));
    Nbp = round(Z/2);
    nn = 1:length(input_signal)+M;


    for i = 1:M-1
        % Method 1
        delay_adc(i) = i/M; % (стр.6,(16)), создаем массив [1:3] из задержанного сигнала ADC0 на различные значения [0.75 0.5 0.25]
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

        %% Algorithm for working in different zones of Nyquist
        nn1 = nn(i:end-(M-i+1)) + delay_adc(i);
        yhil = hilbert(yri);

        cc = cos(2*pi*Nbp*nn1);
        ss = sin(2*pi*Nbp*nn1);

        % plot([real(y_hil(1:100)), imag(y_hil(1:100))]);
        yri_cos = cc.' .* yri; 
        % % plot(yri_cos(1:100))
        yri_sin = ss.' .* yhil; 

        % yri_cos_imag = zeros(length(yri_cos),1);
        % yri_cos_complex = complex(yri_cos,yri_cos_imag);

        if (mod(Z,2) == 0)
            yric = yri_cos + yri_sin;
        else
            yric = yri_cos - yri_sin;
        end

        yri_cut(:,i) = yri((N_taps-1)/2:end); % удаляем переходные процессы (N-1)/2

    end
end