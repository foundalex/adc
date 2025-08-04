function yri_cut = fractional_delays(input_signal, M, N_taps, Z)

    ww = blackman(N_taps);
    Nbp = floor(Z/2);
    nn = 1:length(input_signal(:,1));
    n = (1:1:N_taps);
    del_proc = ((N_taps-1)/2);

    for i = 1:M
        % Method 1
        delay_adc = i/M; % (стр.6,(16)), создаем массив [1:3] из задержанного сигнала ADC0 на различные значения [0.75 0.5 0.25]
        D = ((N_taps-1)/2) + delay_adc;
        hri_m(:,i) = sinc(n - D);

        % plot([ww,hri_m(:,i)])

        % hri_m(:,i) = sin(pi * (n - D)) ./ (pi * (n - D));
        % grpdelay(hri_m(:,i),1,256,1);
        % mean(grpdelay(hri_m(:,i),1,256,1))

        hri_m(:,i) = hri_m(:,i) .* ww;
        hri_m(:,i) = hri_m(:,i) ./ sum(hri_m(:,i)); 

        freqz(hri_m(:,1),1,512, 'whole', 1000000000);  
        yri(:,i) = filter(hri_m(:,i), 1, input_signal(:,1)); % (стр.6 (15))

        % [h,i0,mbw] = designFracDelayFIR(delay_adc,0.5)
        % freqz(h,1,512, 'whole', 1000000000);  

        %% Algorithm for working in different zones of Nyquist
        yhil = hilbert(yri(:,i));
        yhil_imag(:,i) = imag(yhil);
        nn1 = nn + delay_adc;

        cc = cos(2*pi*nn*Nbp);
        ss = sin(2*pi*nn*Nbp);

        yri_cos = cc.' .* yri(:,i); 
        yri_sin = ss.' .* yhil_imag(:,i); 

        if (mod(Z,2) == 0)
            yric(:,i) = yri_cos + yri_sin;
        else
            yric(:,i) = yri_cos - yri_sin;
        end
        yri_cut(:,i) = yric(((N_taps-1)/2)+1:end,i);

    
    end
    % yri_cut(:,1) = input_signal(1:end - ((N_taps-1)/2),1);

    % er = [zeros(36,1); input_signal(1:end-36,2)] - yri(:,1);


    % freqz(hri_m(:,2),1);   
    figure(9);
    plot([[zeros(36,1); input_signal(1:54,2)], yri(1:90,1), yric(1:90,1)]);

    % subplot(8,1,2)
    % plot([[zeros(36,1); input_signal(1:34,3)], yri(1:70,2), yric(1:70,2)]);
    % subplot(8,1,3)
    % plot([[zeros(36,1); input_signal(1:34,4)], yri(1:70,3), yric(1:70,3)]);
    % subplot(8,1,4)
    % plot([[zeros(36,1); input_signal(1:34,5)], yri(1:70,4), yric(1:70,4)]);
    % subplot(8,1,5)
    % plot([[zeros(36,1); input_signal(1:34,6)], yri(1:70,5), yric(1:70,5)]);
    % subplot(8,1,6)
    % plot([[zeros(36,1); input_signal(1:34,7)], yri(1:70,6), yric(1:70,6)]);
    % subplot(8,1,7)
    % plot([[zeros(36,1); input_signal(1:34,8)], yri(1:70,7), yric(1:70,7)]);
    % subplot(2,1,2)
    % plot([zeros(36,1); input_signal(1:54,1)]);


    % spectrumScope = spectrumAnalyzer(SampleRate=1000000000, ...            
    %             AveragingMethod='exponential', ForgettingFactor=0, ...
    %             YLimits=[-30 10],ShowLegend=true, Method='Welch');
    % spectrumScope.WindowLength = 2048;
    % spectrumScope.FrequencyResolutionMethod = "window-length";
    % spectrumScope.PlotAsTwoSidedSpectrum=true;
    % spectrumScope.DistortionMeasurements.Enabled = true;
    % 
    % spectrumScope([input_signal(100:4096+100,2), yric(100:4096+100,1), yric(100:4096+100,2)]);
end