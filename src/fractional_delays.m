function [yri_cut, yri_cut_int, yri_cut1, sig_adc] = fractional_delays(input_signal, input_signal_int, M, N_taps, Z)
    n = (0:1:N_taps-1);
    Nbp = floor(Z/2); % стр 7. (24)
    del_proc = ((N_taps-1)/2);
    nn = 1:length(input_signal(:,1));
    %% Fractional filter coeff
    delay_adc = (1/M:1/M:1); % (стр.6,(16)), создаем массив на различные значения задержек
    w_blackman = 0.42 - 0.5 * cos(2*pi*n/(N_taps-1)) + 0.08 * cos(4*pi*n/(N_taps-1)); % Blackman window

    %% Hilbert filter coeff
    hh = (2./((n-del_proc)*pi)).*(sin(((n-del_proc)*pi)./2)).^2;
    hh(1) = 0;
    hh(37) = 0;
    hh_m = (hh .* w_blackman).';
    hh_m_int = int16(round(hh_m * 2^15))'; % int coeff fi(1,16,15)
    %%


   
    for i = 1:M-1
        w_blackman_fractional(:,i) = 0.42 - 0.5 * cos(2*pi*(n+delay_adc(i))/(N_taps-1)) + 0.08 * cos(4*pi*(n+delay_adc(i))/(N_taps-1)); % shift Blackman window
        D = del_proc - delay_adc(i); % delay (N-1)/2 - d = causal filter
        hri_m(:,i) = sinc(n-D); % shift impulse response on D = Dint - d for fractional delay filter
        %% coeff fractional filter
        hri_m(:,i) = hri_m(:,i) .* w_blackman_fractional(:,i); 
        % hri_m(:,i) = hri_m(:,i) ./ sum(hri_m(:,i));
        hrim_fi(:,i) = fi(hri_m(:,i),1,19,18);
        hrim_int(:,i) = int32(round(hrim_fi(:,i) * 2^18)); % fi(1,19,18)

        %% filter
        yri = filter(hri_m(:,i), 1, input_signal(:,1)); % (стр.6 (15))

        %% codegen
        % buildInstrumentedMex  fir_filter -coder -o fir_filter -args {hrim_int, input_signal_int(:,1)} -histogram
    
        y_out_arrayInt = fir_filter(hrim_int(:,i), input_signal_int(:,1)); % (стр.6 (15)) fi(1,19,18) * fi(1,12,11) = fi(1,31,29)
        y_out_array_int16 = int16(round(y_out_arrayInt/262144)); % fi(1,31,29) - 19 = fi(1,12,10)

        % showInstrumentationResults fir_filter

        % clear InstrumentationResults fir_filter
        % clear fir_filter
        % delete fir_filter.mexw64
        % delete codegen/mex/fir_filter/fir_filter.mexw64

        % [yy(:,i), ff] = freqz(hri_m(:,i),1,1024, 'whole', 1000000000);
        % y_out_array_int = int32(filter(hrim_int(:,i), 1, input_signal_int(:,1))); % (стр.6 (15)) fi(1,19,18) * fi(1,12,11) = fi(1,31,29)

        %% Algorithm for working in different zones of Nyquist
        %%

        ymi = filter(hh_m.', 1, yri);
        % ymi_int = int32((filter(hh_m_int, 1, y_out_array_int16(:,i)))); % (стр.6 (15)) % fi(1,16,15) * fi(1,14,12) = fi(1,30,27)
        ymi_arrayInt = fir_filter(int32(hh_m_int), y_out_array_int16); % (стр.6 (15)) fi(1,16,15) * fi(1,14,11) = fi(1,30,27)
        ymi_int16 = int16(round(ymi_arrayInt/32768)); % fi(1,30,27) - 16 bit = fi(1,14,11)

        % figure(2)
        % plot([ymi(1:500), double(ymi_int16(1:500))*2^-11]);

        %%
        yhil_imag = [ymi(del_proc+1:end); zeros(del_proc,1)];
        yhil_imag_int = [ymi_int16(del_proc+1:end); zeros(del_proc,1)];

        nn1 = nn + delay_adc(i);
        a1 = 2*pi*nn1*Nbp;
        cc = cos(a1);
        ss = sin(a1);

        %%
        if Z == 2 | Z == 3
           if (delay_adc(i) == 0.25)
                yri_cos_int16 = int16(0);
                yri_sin_int16 = yhil_imag_int;
           elseif (delay_adc(i) == 0.5)
                yri_cos_int16 = -y_out_array_int16;
                yri_sin_int16 = int16(0);
           elseif (delay_adc(i) == 0.75)
                yri_cos_int16 = int16(0);
                yri_sin_int16 = -yhil_imag_int;
           end
        elseif Z == 4
           if (delay_adc(i) == 0.25)
                yri_cos_int16 = -y_out_array_int16;
                yri_sin_int16 = int16(0);
           elseif (delay_adc(i) == 0.5)
                yri_cos_int16 = y_out_array_int16;
                yri_sin_int16 = int16(0);
           elseif (delay_adc(i) == 0.75)
                yri_cos_int16 = -y_out_array_int16;
                yri_sin_int16 = int16(0);
           end
        else 
            yri_cos_int16 = y_out_array_int16;
            yri_sin_int16 = int16(0);
        end


        yri_cos = cc' .* yri; 
        yri_sin = ss' .* yhil_imag; 

        % figure(4);
        % plot([yri_cos(1:200), double(yri_cos_int16(1:200)), yri_sin(1:200), double(yri_sin_int16(1:200))]);

        %%
        if (mod(Z,2) == 0)
            yric(:,i) = yri_cos + yri_sin;
            yric_int16(:,i) = yri_cos_int16 + yri_sin_int16; % fi(1,12,10) + fi(1,14,11) = fi(1,15,11) 

        else
            yric(:,i) = yri_cos - yri_sin;
            yric_int16(:,i) = yri_cos_int16 - yri_sin_int16; % fi(1,12,10) - fi(1,14,11) = fi(1,15,11) 
        end

        % figure(4);
        % plot([yric(1:200,i), double(yric_int16(1:200,i))*2^-11]);

        yri_cut(:,i+1) = yric(del_proc+1:end,i);
        yri_cut_int(:,i+1) = yric_int16(del_proc+1:end,i);
    
    end
    
    yri_cut(:,1) = input_signal(1:end-del_proc,1);
    yri_cut_int(:,1) = input_signal_int(1:end-del_proc,1);


    yri_cut(end-del_proc:end,:) = [];
    yri_cut_int(end-del_proc:end,:) = [];

    % figure(9);
    % % subplot(2,1,1)
    % % plot([[zeros(del_proc,1); input_signal(1:274,2)], yric(1:del_proc+274,1)]);
    % % subplot(3,1,2)
    % plot([input_signal(1:274,1), yri_cut(1:274,1)]);
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

    %% test signal after fractional delay filters
    for j = 1:M
        yri_cut1(:,j) = double(yri_cut_int(:,j))*2^-11;
    end

	sig_adc = zeros(M*length(yri_cut(:,1)),1);
    sig_adc_int = zeros(M*length(yri_cut_int(:,1)),1);
    
	for i = 1:M
		sig_adc(i:M:end) = yri_cut(:,i);
        sig_adc_int(i:M:end) = yri_cut1(:,i);
    end

    % figure(4);
    % plot([sig_adc(1:250), sig_adc_int(1:250)]);
    % % 
    % figure(10);
    % subplot(2,1,1);
    % snr(sig_adc, 8000000000);
    % subplot(2,1,2);
    % snr(sig_adc_int, 8000000000);

end