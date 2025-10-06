function [yri_cut, yri_cut_int, yri_cut1, sig_adc] = fractional_delays(input_signal, input_signal_int, M, N_taps, Z)
    n = (0:1:N_taps-1);
    Nbp = floor(Z/2); % стр 7. (24)
    del_proc = ((N_taps-1)/2);
    nn = 1:length(input_signal(:,1));
    w_blackman = 0.42 - 0.5 * cos(2*pi*n/(N_taps-1)) + 0.08 * cos(4*pi*n/(N_taps-1)); % Blackman window
    %% Fractional filter coeff
    %%
    delay_adc = (1/M:1/M:1); % (стр.6,(16)), создаем массив на различные значения задержек
    D = del_proc - delay_adc; % delay (N-1)/2 - d = causal filter
    hri_m = sinc(n'- D); % shift impulse response on D = Dint - d for fractional delay filter
    w_blackman_fractional = 0.42 - 0.5 * cos(2*pi*(n'+delay_adc)/(N_taps-1)) + 0.08 * cos(4*pi*(n'+delay_adc)/(N_taps-1)); % shift Blackman window
    hri_w = hri_m .* w_blackman_fractional; 

    %% integer
    fractional_width = 13;
    hrim_fi = fi(hri_w, 1,fractional_width+1,fractional_width);
    hrim_int = int16(hrim_fi * 2^fractional_width); % fi(1,6,5)

    % figure(2);
    % freqz(double(hrim_int)*2^-fractional_width,1,1024, 'whole', 1000000000);


    %% Hilbert filter coeff
    %%
    hh = (2./((n-del_proc)*pi)).*(sin(((n-del_proc)*pi)./2)).^2;
    hh(1) = 0;
    hh(37) = 0;
    hh_m = (hh .* w_blackman).';


    hilbert_width = 13;
    

    hh_m_fi = fi(hh_m,1,hilbert_width,hilbert_width-1);
    hh_m_int = int16(round(hh_m_fi * 2^(hilbert_width-1))); % int coeff fi(1,4,3)

    % [y, f] = freqz(hh_m, 1, 1024, 'whole', 1000000000);
    % [y1,f1] = freqz(double(hh_m_int)*2^-(hilbert_width-1), 1, 1024, 'whole', 1000000000);
    % figure(2);
    % plot(f,abs(y), f1, abs(y1));

    %%


   
    for i = 1:M-1
        




        yri = filter(hri_w(:,i), 1, input_signal(:,1)); % filter (стр.6 (15))

        %% codegen
        % buildInstrumentedMex  fir_filter -coder -o fir_filter -args {hrim_int, input_signal_int(:,1)} -histogram
    
        y_fractional_outInt = fir_filter(hrim_int(:,i), input_signal_int(:,1)); % (стр.6 (15)) 
        y_fractional_outInt16 = int16(bitshift(y_fractional_outInt, -width_fractional)); % fi(1,13,11)

        figure(4);
        plot([yri(1:250), double(y_fractional_outInt16(1:250))*2^0]);
        % 
        figure(10);
        subplot(2,1,1);
        snr(yri, 8000000000);
        subplot(2,1,2);
        snr(double(y_fractional_outInt16)*2^0, 8000000000);


        y_fractional_outInt_max(:,i) = max(y_fractional_outInt16);
        y_fractional_outInt_min(:,i) = min(y_fractional_outInt16);

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
        ymi_HilbertInt = fir_filter(hh_m_int, y_fractional_outInt16); % (стр.6 (15)) fi(1,16,15) * fi(1,14,11) = fi(1,30,27)
        % ymi_int16 = int16(round(ymi_HilbertInt/32768)); % fi(1,30,27) - 16 bit = fi(1,14,11)
        ymi_int16 = ymi_HilbertInt; % (bitshift(ymi_HilbertInt, -1)); % fi(1,13,11)
        % figure(2)
        % plot([ymi(1:500), double(ymi_int16(1:500))*2^-11]);

        figure(4);
        plot([ymi(1:250), double(ymi_HilbertInt(1:250))*2^-12]);
        % 
        figure(10);
        subplot(2,1,1);
        snr(ymi, 8000000000);
        subplot(2,1,2);
        snr(double(ymi_HilbertInt)*2^-12, 8000000000);

        y_out_Hilbert_max(:,i) = max(ymi_HilbertInt);
        y_out_Hilbert_min(:,i) = min(ymi_HilbertInt);

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
                yri_cos_int16 = -y_fractional_outInt16;
                yri_sin_int16 = int16(0);
           elseif (delay_adc(i) == 0.75)
                yri_cos_int16 = int16(0);
                yri_sin_int16 = -yhil_imag_int;
           end
        elseif Z == 4
           if (delay_adc(i) == 0.25)
                yri_cos_int16 = -y_fractional_outInt16;
                yri_sin_int16 = int16(0);
           elseif (delay_adc(i) == 0.5)
                yri_cos_int16 = y_fractional_outInt16;
                yri_sin_int16 = int16(0);
           elseif (delay_adc(i) == 0.75)
                yri_cos_int16 = -y_fractional_outInt16;
                yri_sin_int16 = int16(0);
           end
        else 
            yri_cos_int16 = y_fractional_outInt16;
            yri_sin_int16 = int16(0);
        end

        yri_cos = cc' .* yri; 
        yri_sin = ss' .* yhil_imag; 

        figure(5);
        plot([yri_cos(1:200), double(yri_cos_int16(1:200))]);

        figure(10);
        subplot(2,1,1);
        snr(yri_cos, 8000000000);
        subplot(2,1,2);
        snr(double(yri_cos_int16), 8000000000);


        %%
        if (mod(Z,2) == 0)
            yric(:,i) = yri_cos + yri_sin;
            yric_int16(:,i) = yri_cos_int16 + yri_sin_int16; % fi(1,12,10) + fi(1,14,11) = fi(1,15,11) 

        else
            yric(:,i) = yri_cos - yri_sin;
            yric_int16(:,i) = yri_cos_int16 - yri_sin_int16; % fi(1,12,10) - fi(1,14,11) = fi(1,15,11) 
        end

        % figure(4);
        % plot([yric(1:200,i), double(yric_int16(1:200,i))*2^-14]);
        % 
        % figure(10);
        % subplot(2,1,1);
        % snr(yric, 8000000000);
        % subplot(2,1,2);
        % snr(double(yric_int16)*2^-14, 8000000000);


        yri_cut(:,i+1) = yric(del_proc+1:end,i);
        yri_cut_int(:,i+1) = yric_int16(del_proc+1:end,i);
    
    end
    
    yri_cut(:,1) = input_signal(1:end-del_proc,1);
    yri_cut_int(:,1) = input_signal_int(1:end-del_proc,1);

    yri_cut(end-del_proc:end,:) = [];
    yri_cut_int(end-del_proc:end,:) = [];

    % figure(9);
    % subplot(4,1,1)
    % plot([yri_cut(1:274,1), double(yri_cut_int(1:274,1))*2^-11]);
    % subplot(4,1,2)
    % plot([yri_cut(1:274,2), double(yri_cut_int(1:274,2))*2^-16]);
    % subplot(4,1,3)
    % plot([yri_cut(1:274,3), double(yri_cut_int(1:274,3))*2^-16]);
    % subplot(4,1,4)
    % plot([yri_cut(1:274,4), double(yri_cut_int(1:274,4))*2^-16]);


    %% test signal after fractional delay filters

	sig_adc = zeros(M*length(yri_cut(:,1)),1);
    sig_adc_int = zeros(M*length(yri_cut_int(:,1)),1);
    
    % yri_cut1(:,1) = double(yri_cut_int(:,1))*2^-16;
    
    % yri_cut1(:,1) = double(yri_cut_int(:,1))*2^-11;
    % yri_cut1(:,2) = double(yri_cut_int(:,2))*2^-16;
    % yri_cut1(:,3) = double(yri_cut_int(:,3))*2^-16;
    % yri_cut1(:,4) = double(yri_cut_int(:,4))*2^-16;

	for i = 1:M
		sig_adc(i:M:end) = yri_cut(:,i);
        sig_adc_int(i:M:end) = double(yri_cut_int(:,i));
    end

    figure(4);
    plot([sig_adc(1:250), sig_adc_int(1:250)]);
    % 
    figure(10);
    subplot(2,1,1);
    snr(sig_adc, 8000000000);
    subplot(2,1,2);
    snr(sig_adc_int, 8000000000);

    max_fractional = max(y_fractional_outInt_max);
    min_fractional = min(y_fractional_outInt_min);

end