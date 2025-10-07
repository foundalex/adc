
% 1) Hu.M, Yi.P, (2022), Digital Calibration for Gain, Time Skew, and Bandwidth Mismatch 
%    in Under-Sampling Time-Interleaved System
% 2) Джиган В.И, Адаптивные фильтры
% 3) Айфичер Э, Джервис Б, Цифровая обработка сигналов. Практический подход
% 4) Behrouz Farhang-Boroujeny, Adaptive Filters Theory and Applications 

function [x_after_adc, x_after_adc_int, snr_s] = adc_calibration(sim_options, adc_input, s_to_subadc, s_after_subadc, Z)

    %% Calibration algorithm 1 (Fractional delays)
    n = (0:1:sim_options.N-1);
    Nbp = floor(Z/2); % стр 7. (24)
    del_proc = ((sim_options.N-1)/2);
    nn = 1:length(adc_input(:,1));
    w_blackman = 0.42 - 0.5 * cos(2*pi*n/(sim_options.N-1)) + 0.08 * cos(4*pi*n/(sim_options.N-1)); % Blackman window
    %% Fractional filter coeff
    %%
    delay_adc = (1/sim_options.M:1/sim_options.M:1); % (стр.6,(16)), создаем массив на различные значения задержек
    D = del_proc - delay_adc; % delay (N-1)/2 - d = causal filter
    hri_m = sinc(n'- D); % shift impulse response on D = Dint - d for fractional delay filter
    w_blackman_fractional = 0.42 - 0.5 * cos(2*pi*(n'+ delay_adc)/(sim_options.N-1)) + 0.08 * cos(4*pi*(n'+ delay_adc)/(sim_options.N-1)); % shift Blackman window
    hri_w = hri_m .* w_blackman_fractional; 

    fractional_width = [16, 19, 21];

    [y, f] = freqz(hri_w(:,1), 1,1024, 'whole', 1000000000);
    for k = 1:sim_options.M-1
        hrim_fi_test = fi(hri_w(:,1), 1,fractional_width(k),fractional_width(k)-1);
        hrim_int_test = int32(hrim_fi_test * 2^(fractional_width(k)-1));
        [y1(:,k), f1(:,k)] = freqz(double(hrim_int_test)*2^-(fractional_width(k)-1),1,1024, 'whole', 1000000000);
    end

    % figure(2);
    % plot(f, abs(y), f, abs(y1(:,1)), f, abs(y1(:,2)), f, abs(y1(:,3)));
    % title('Влияние разрядностей коэффициентов на АЧХ фильтра дробной задержки')
    % xlabel('Частота') 
    % ylabel('Коэффициент передачи') 
    % legend({'double','16 бит', '19 бит', '21 бит'},'Location','northeast')

    %% Hilbert filter coeff
    %%
    hh = (2./((n-del_proc)*pi)).*(sin(((n-del_proc)*pi)./2)).^2;
    hh(1) = 0;
    hh(37) = 0;
    hh_m = (hh .* w_blackman).';

    hilbert_width = 13;
    


    nn1 = nn' + delay_adc;
    a1 = 2*pi*nn1*Nbp;
    c_os = cos(a1);
    s_os = sin(a1);

    for r = 1:1
	    % Fractional delays of ADC0 signal
        for i = 1:sim_options.M-1
	        [yri, ymi, y_fractional_outInt, ymi_HilbertInt] = fractional_delays(adc_input(:,1), adc_input(:,1), Z, hri_w(:,i), hh_m, fractional_width(2), hilbert_width);

            max_fractional(:,i) = max(y_fractional_outInt);
            min_fractional(:,i) = min(y_fractional_outInt);

            yhil_imag = [ymi(del_proc+1:end); zeros(del_proc,1)];
            yhil_imag_int = [ymi_HilbertInt(del_proc+1:end); zeros(del_proc,1)];

            %% Algorithm for working in different zones of Nyquist
            yri_cos = c_os(:,i) .* yri; 
            yri_sin = s_os(:,i) .* yhil_imag; 

            if Z == 2 | Z == 3
                if (delay_adc(i) == 0.25)
                    yri_cos_int16 = int32(0);
                    yri_sin_int16 = yhil_imag_int(:,i);
                elseif (delay_adc(i) == 0.5)
                    yri_cos_int16 = -y_fractional_outInt;
                    yri_sin_int16 = int32(0);
                elseif (delay_adc(i) == 0.75)
                    yri_cos_int16 = int32(0);
                    yri_sin_int16 = -yhil_imag_int(:,i);
                end
            elseif Z == 4
                if (delay_adc(i) == 0.25)
                    yri_cos_int16 = -y_fractional_outInt;
                    yri_sin_int16 = int32(0);
                elseif (delay_adc(i) == 0.5)
                    yri_cos_int16 = y_fractional_outInt;
                    yri_sin_int16 = int32(0);
                elseif (delay_adc(i) == 0.75)
                    yri_cos_int16 = -y_fractional_outInt;
                    yri_sin_int16 = int32(0);
                end
            else 
                yri_cos_int16 = y_fractional_outInt;
                yri_sin_int16 = int32(0);
            end


            if (mod(Z,2) == 0)
                yric = yri_cos + yri_sin;
                yric_int16 = yri_cos_int16 + yri_sin_int16; % fi(1,12,10) + fi(1,14,11) = fi(1,15,11) 
            else
                yric = yri_cos - yri_sin;
                yric_int16 = yri_cos_int16 - yri_sin_int16; % fi(1,12,10) - fi(1,14,11) = fi(1,15,11) 
            end
            %%

            yri_cut(:,i+1) = yric(del_proc+1:end);
            yri_cut_int(:,i+1) = yric_int16(del_proc+1:end);

        end


        yri_cut(:,1) = adc_input(1:end-del_proc,1);
        yri_cut_int(:,1) = adc_input(1:end-del_proc,1);

        yri_cut(end-del_proc:end,:) = [];
        yri_cut_int(end-del_proc:end,:) = [];


        %% test signal after fractional delay filters
        sig_adc = zeros(sim_options.M*length(yri_cut(:,1)),1);
        sig_adc_int = zeros(sim_options.M*length(yri_cut_int(:,1)),1);


        yri_cut1(:,1) = double(yri_cut_int(:,1));
        yri_cut1(:,2) = double(yri_cut_int(:,2))*2^-(fractional_width(2)-1);
        yri_cut1(:,3) = double(yri_cut_int(:,3))*2^-(fractional_width(2)-1);
        yri_cut1(:,4) = double(yri_cut_int(:,4))*2^-(fractional_width(2)-1);


	    for i = 1:sim_options.M
            sig_adc(i:sim_options.M:end) = yri_cut(:,i);
            sig_adc_int(i:sim_options.M:end) = double(yri_cut1(:,i));
        end

        % figure(4);
        % plot([sig_adc(1:250), sig_adc_int(1:250)]);

        snr_s(r+1) = snr(sig_adc_int, sim_options.Fs/sim_options.Inter);

        figure(10);
        subplot(2,1,1);
        snr(sig_adc, 8000000000);
        subplot(2,1,2);
        snr(sig_adc_int, 8000000000);

        % yri_cut = [];
        % yri_cut_int = [];

    end


    snr_s(1) = snr(sig_adc, sim_options.Fs/sim_options.Inter);


    % figure(3);
    % plot(sim_options.freq, snr_s(1), '-o', sim_options.freq, snr_s(2), '-o', sim_options.freq, snr_s(3), '-o', sim_options.freq, snr_s(4), '-o');
    % title('SNR (dB)')
    % xlabel('Частота') 
    % ylabel('SNR (dB)') 
    % legend({'double', '16 бит', '19 бит', '21 бит'}, 'Location','northwest');





    max_l = max(max_fractional);
    min_l = min(min_fractional);


    x_after_adc = 0;
    x_after_adc_int = 0;

    %% Calibration algorithm 2 (Least Mean Squares)

    adc_input = double(adc_input)*2^-11;
    yri_cut = double(yri_cut_int)*2^-11; 

    % adc_input = double(adc_input_int);
    % yri_cut = double(yri_cut_int);

    [y_array, y_array_int] = least_mean_squares(adc_input, adc_input, yri_cut, yri_cut_int, sim_options.M, sim_options.N1, sim_options.Width);

    % create main signal after LS algorithm (switch after sub-adc)
    x_after_adc = zeros(length(y_array)*sim_options.M,1);
    x_after_adc_int = zeros(length(y_array_int)*sim_options.M,1);

    for i = 1:sim_options.M
        if i == 1
            x_after_adc(i:sim_options.M:end) = yri_cut1(1:length(y_array),1);
            % x_after_adc(i:sim_options.M:end) = double(yri_cut_int(1:length(y_array),1));
            x_after_adc_int(i:sim_options.M:end) = double(yri_cut_int(1:length(y_array),1));
        else
            x_after_adc(i:sim_options.M:end) = y_array(:,i-1);
            x_after_adc_int(i:sim_options.M:end) = double(y_array_int(:,i-1)) * 2^-sim_options.Width;
        end
    end
    % 
    % figure(3);
    % subplot(2,1,1)
    % plot([x_after_adc_int]);
    % title('Выход адаптивного фильтра int')
    % xlabel('Номер отсчета') 
    % ylabel('Амплитуда') 
    % 
    % subplot(2,1,2)
    % plot([x_after_adc]);
    % title('Выход адаптивного фильтра double')
    % xlabel('Номер отсчета') 
    % ylabel('Амплитуда');






    % legend('double', 'double', 'Исходный сигнал с ошибками')

    % figure(4);
    % subplot(4,1,1);
    % snr(s_to_subadc(1:length(s_to_subadc)), sim_options.Fs/sim_options.Inter);
    % subplot(4,1,2);
    % snr(s_after_subadc(1:length(s_after_subadc)), sim_options.Fs/sim_options.Inter);
    % subplot(4,1,3);
    % snr(x_after_adc(1:length(x_after_adc)), sim_options.Fs/sim_options.Inter);
    % subplot(4,1,4);
    % snr(x_after_adc_int(1:length(x_after_adc_int)), sim_options.Fs/sim_options.Inter);

end
