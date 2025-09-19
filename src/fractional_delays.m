function [yri_cut, yri_cut_int, yri_cut1, sig_adc] = fractional_delays(input_signal, input_signal_int, M, N_taps, Z)

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
        % hri_m(:,i) = hri_m(:,i) ./ sum(hri_m(:,i));

        %% integer
        hrim_fi(:,i) = fi(hri_m(:,i),1,19,18);
        hrim_int(:,i) = int32(round(hrim_fi(:,i) * 2^18)); % fi(1,19,18)
        % figure(4);
        % plot([hri_m(:,i), double(hrim_fi(:,i)), double(hrim_int(:,i)) * 2^-18]);

        %%
        [yy(:,i), ff] = freqz(hri_m(:,i),1,1024, 'whole', 1000000000);
        yri(:,i) = filter(hri_m(:,i), 1, input_signal(:,1)); % (стр.6 (15))
        %% integer
        filt_width = 12;
        y_out_array_int(:,i) = int32(filter(hrim_int(:,i), 1, input_signal_int(:,1))); % (стр.6 (15)) fi(1,19,18) * fi(1,12,11) = fi(1,31,29)
        y_out_array_int16(:,i) = (int16(round(y_out_array_int(:,i)/131072))).'; % fi(1,31,29) - 17 = fi(1,14,12)

        % figure(4);
        % plot([yri(1:200,i), double(y_out_array_int(1:200,i)) * 2^-29, double(y_out_array_int16(1:200,i)) * 2^-filt_width]);

        % figure(3);
        % subplot(3,1,1);
        % sfdr(yri(:,i), 1000000000);
        % subplot(3,1,2);
        % sfdr(double(y_out_array_int(:,i)) * 2^-28, 1000000000);
        % subplot(3,1,3);
        % sfdr(double(y_out_array_int16(:,i)) * 2^-filt_width, 1000000000);

        %% Algorithm for working in different zones of Nyquist

        hh = (2./((n-del_proc)*pi)).*(sin(((n-del_proc)*pi)./2)).^2;
        hh(1) = 0;
        hh(37) = 0;

        hh_m = (hh .* w_blackman).';
        %% integer

        hh_m_int = int16(round(hh_m * 2^15)); % fi(1,16,15)

        filt_width_1 = 11;
        ymi = filter(hh_m.', 1, yri(:,i));
        ymi_int = int32((filter(hh_m_int, 1, y_out_array_int16(:,i)))); % (стр.6 (15)) % fi(1,16,15) * fi(1,14,12) = fi(1,30,27)

        ymi_int16 = int16(round(ymi_int/65536)); % fi(1,30,27) - 16 bit = fi(1,14,11)

        % figure(4);
        % plot([ymi(1:200), double(ymi_int(1:200)) * 2^-27, double(ymi_int16(1:200)) * 2^-filt_width_1]);
        % 
        % figure(3);
        % subplot(3,1,1);
        % snr(ymi, 1000000000);
        % subplot(3,1,2);
        % snr(double(ymi_int) * 2^-27, 1000000000);
        % subplot(3,1,3);
        % snr(double(ymi_int16) * 2^-filt_width_1, 1000000000);

        %%
        yhil_imag(:,i) = [ymi(del_proc+1:end); zeros(del_proc,1)];
        yhil_imag_int(:,i) = [ymi_int16(del_proc+1:end); zeros(del_proc,1)];

        nn1 = nn + delay_adc(i);
        cc = cos(2*pi*nn1*Nbp);
        ss = sin(2*pi*nn1*Nbp);

        %%
        cc_fix = fi(cc,1,16,15);
        ss_fix = fi(ss,1,16,15);

        % figure(4);
        % plot([cc(1:200), double(cc_fix(1:200)), ss(1:200), double(ss_fix(1:200))]);
        %%
        cc_int = int16(round(cc_fix * 2^15));
        ss_int = int16(round(ss_fix * 2^15));
        %%
        yri_cos = cc.' .* yri(:,i); 
        yri_sin = ss.' .* yhil_imag(:,i); 
        %% 
        yri_cos_int = int32(cc_int.') .* int32(y_out_array_int16(:,i)); % fi(1,16,15) * fi(1,14,12) = fi(1,30,27);
        yri_cos_int16 = int16(round(yri_cos_int/65536)); % fi(1,31,27) - 16 bit = fi(1,15,11)

        % figure(4);
        % plot([yri_cos(1:200), double(yri_cos_int(1:200)) * 2^-27, double(yri_cos_int16(1:200))*2^-11]);

        yri_sin_int = int32(ss_int.') .* int32(yhil_imag_int(:,i)); % fi(1,16,15) + fi(1,14,11) = fi(1,30,26);
        yri_sin_int16 = int16(round(yri_sin_int/32768)); % fi(1,30,26) - 15 bit = fi(1,15,11)

        % figure(4);
        % plot([yri_sin(1:200), double(yri_sin_int16(1:200))*2^-10]);

        %%
        if (mod(Z,2) == 0)
            yric(:,i) = yri_cos + yri_sin;
            yric_int16(:,i) = int32(yri_cos_int16) + int32(yri_sin_int16); % fi(1,15,11) + fi(1,15,11) = fi(1,16,11) 

            yric_double = double(yric_int16)*2^-11;

            % figure(4);
            % plot([yric(1:200), yric_double(1:200)]);

        else
            yric(:,i) = yri_cos - yri_sin;
            yric_int16(:,i) = yri_cos_int16 - yri_sin_int16; % fi(1,14,10) - fi(1,14,10) 

            yric_double = double(yric_int16)*2^-11;
            % 
            % figure(4);
            % plot([yric(1:200), yric_double(1:200)]);
        end

        yri_cut(:,i+1) = yric(del_proc+1:end,i);
        yri_cut_int(:,i+1) = yric_int16(del_proc+1:end,i);
    
    end
    yri_cut(:,1) = input_signal(1:end-del_proc,1);
    yri_cut_int(:,1) = input_signal_int(1:end-del_proc,1);


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

    %% test signal after fractional delay filters
    for j = 1:M
        yri_cut1(:,j) = double(yri_cut_int(:,j))*2^-11;
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
    plot([sig_adc, sig_adc_int]);
    % 
    figure(10);
    subplot(2,1,1);
    snr(sig_adc, 8000000000);
    subplot(2,1,2);
    snr(sig_adc_int, 8000000000);
    


end