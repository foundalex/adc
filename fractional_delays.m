function yri_cut = fractional_delays(input_signal, M, N_taps, Z)

    n = (0:1:N_taps-1);
    Nbp = floor(Z/2); % стр 7. (24)
    nn = 1:length(input_signal(:,1));
    del_proc = ((N_taps-1)/2);
    delay_adc = (1/M:1/M:1); % (стр.6,(16)), создаем массив из задержанного сигнала ADC0 на различные значения

    for i = 1:M-1
        % Method 1
        D = del_proc - delay_adc(i); % delay (N-1)/2 - d = causal filter

        w_blackman_fractional = 0.42 - 0.5 * cos(2*pi*(n+delay_adc(i))/(N_taps-1)) + 0.08 * cos(4*pi*(n+delay_adc(i))/(N_taps-1)); % shift Blackman window
        w_blackman = 0.42 - 0.5 * cos(2*pi*n/(N_taps-1)) + 0.08 * cos(4*pi*n/(N_taps-1)); % Blackman window

        % hri_m(:,i) = sinc(n-D) - 0.1*sinc(0.1*(n-D)); % try high pass filter

        hri_m(:,i) = sinc(n-D); % shift impulse response on D = Dint - d for fractional delay filter

        %% for example
        hri_m1(:,i) = sinc(n-del_proc); % shift impulse response on Dint
       
        %%
        % plot(n.',w1.', n.',w.', n.',hri_m(:,i), n.',hri_m1(:,i))
        % grpdelay(hri_m(:,i),1,256,1);
        % mean(grpdelay(hri_m(:,i),1,256,1))

        hri_m(:,i) = hri_m(:,i) .* w_blackman_fractional.'; 
        % hri_m(:,i) = hri_m(:,i) ./ sum(hri_m(:,i)); 

        % freqz(hri_m(:,i),1,1024, 'whole', 1000000000)  
     
        yri(:,i) = filter(hri_m(:,i), 1, input_signal(:,1)); % (стр.6 (15))

        % hri_m(:,i) = firlp2hp(hri_m(:,i)); 
        % freqz(hri_m(:,i),1,1024, 'whole', 1000000000)  
        %%
        % Method 2
        % [h,i0,mbw] = designFracDelayFIR(delay_adc(i),0.98);
        % [hh1] = freqz(h,1,512, 'whole', 1000000000); 
        % yri(:,i) = filter(h, 1, input_signal(:,1)); % (стр.6 (15))
        % plot(ww, abs(hh), ww, abs(hh1))

        %% Algorithm for working in different zones of Nyquist
        % yhil = hilbert(yri(:,i)); 
        % yhil_imag(:,i) = imag(yhil);

        hh = (2./((n-del_proc)*pi)).*(sin(((n-del_proc)*pi)./2)).^2;
        hh(1)= 0;
        hh(37)= 0;
        % plot(hh)
        hh_m = hh .* w_blackman;
        ymi = filter(hh_m.', 1, yri(:,i));

        % plot([yri(1:100,i), ymi(37:100+36), yhil_imag(1:100)]);

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
  
    % figure(9);
    % subplot(9,1,1)
    % plot([[zeros(del_proc,1); input_signal(1:274,2)], yric(1:del_proc+274,1)]);
    % subplot(9,1,2)
    % plot([input_signal(1:274,1), yri_cut(1:274,1)]);
    % subplot(9,1,3)
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