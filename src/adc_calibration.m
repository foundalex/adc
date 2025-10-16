
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

    % [y, f] = freqz(hri_w(:,1), 1,1024, 'whole', 1000000000);
    % for k = 1:sim_options.M-1
    %     hrim_fi_test = fi(hri_w(:,k), 1,fractional_width(1),fractional_width(1)-1);
    %     hrim_int_test = int32(hrim_fi_test * 2^(fractional_width(1)-1));
    %     [y1(:,k), f1(:,k)] = freqz(double(hrim_int_test)*2^-(fractional_width(1)-1),1,1024, 'whole', 1000000000);
    % end
    % 
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
    %% zones Nyquist
    nn1 = nn' + delay_adc;
    a1 = 2*pi*nn1*Nbp;
    c_os = cos(a1);
    s_os = sin(a1);

    for r = 1:1
	    % Fractional delays of ADC0 signal
        for i = 1:sim_options.M-1
	        [yri, ymi, y_fractional_outInt, ymi_HilbertInt] = fractional_delays(adc_input(:,1), Z, hri_w(:,i), hh_m, fractional_width(2), hilbert_width);

            % y_fractional_outInt = fi(1,33,18)                            ymi_HilbertInt = fi(1,62,30)

            max_fractional(:,i) = max(y_fractional_outInt);
            min_fractional(:,i) = min(y_fractional_outInt);

            max_hilbert(:,i) = max(ymi_HilbertInt);
            min_hilbert(:,i) = min(ymi_HilbertInt);

            yhil_imag = [ymi(del_proc+1:end); zeros(del_proc,1)];
            yhil_imag_int = [ymi_HilbertInt(del_proc+1:end); zeros(del_proc,1)];

            %% Algorithm for working in different zones of Nyquist
            yri_cos = c_os(:,i) .* yri; 
            yri_sin = s_os(:,i) .* yhil_imag; 

            if Z == 2 | Z == 3
                if (delay_adc(i) == 0.25)
                    yric_int = yhil_imag_int;
                elseif (delay_adc(i) == 0.5)
                    yric_int = -y_fractional_outInt;
                elseif (delay_adc(i) == 0.75)
                    yric_int = -yhil_imag_int;
                end
            elseif Z == 4
                if (delay_adc(i) == 0.25)
                    yric_int = -y_fractional_outInt;
                elseif (delay_adc(i) == 0.5)
                    yric_int = y_fractional_outInt;
                elseif (delay_adc(i) == 0.75)
                    yric_int = -y_fractional_outInt;
                end
            else 
                yric_int = y_fractional_outInt;
            end


            if (mod(Z,2) == 0)
                yric = yri_cos + yri_sin;
                % yric_int = yri_cos_int + yri_sin_int; % fi(1,63,0)
            else
                yric = yri_cos - yri_sin;
                % yric_int = yri_cos_int - yri_sin_int; % fi(1,63,0)
            end
            %%

            yri_cut(:,i+1) = yric(del_proc+1:end);
            yri_cut_int(:,i+1) = yric_int(del_proc+1:end);

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

        figure(4);
        plot([sig_adc(1:250), sig_adc_int(1:250)]);
        figure(10);
        subplot(2,1,1);
        snr(sig_adc, 8000000000);
        subplot(2,1,2);
        snr(sig_adc_int, 8000000000);

        snr_s(r+1) = snr(sig_adc_int, sim_options.Fs/sim_options.Inter);

        % yri_cut = [];
        % yri_cut_int = [];

    end

    snr_s(1) = snr(sig_adc, sim_options.Fs/sim_options.Inter);

    max_l = max(max_fractional);
    min_l = min(min_fractional);

    max_h = max(max_hilbert);
    min_h = min(min_hilbert);

    %% Calibration algorithm 2 (Least Mean Squares)

    [y_array, y_array_int] = least_mean_squares(double(adc_input), adc_input, yri_cut, yri_cut_int, sim_options.M, sim_options.N1, sim_options.Width);

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
            x_after_adc_int(i:sim_options.M:end) = double(y_array_int(:,i-1)) * 2^-32;
        end
    end
    % 
    figure(3);
    subplot(2,1,1)
    plot([x_after_adc_int]);
    title('Выход адаптивного фильтра int')
    xlabel('Номер отсчета') 
    ylabel('Амплитуда') 

    subplot(2,1,2)
    plot([x_after_adc]);
    title('Выход адаптивного фильтра double')
    xlabel('Номер отсчета') 
    ylabel('Амплитуда');


    % legend('double', 'double', 'Исходный сигнал с ошибками')

    figure(4);
    subplot(4,1,1);
    snr(s_to_subadc(1:length(s_to_subadc)), sim_options.Fs/sim_options.Inter);
    subplot(4,1,2);
    snr(s_after_subadc(1:length(s_after_subadc)), sim_options.Fs/sim_options.Inter);
    subplot(4,1,3);
    snr(x_after_adc(1:length(x_after_adc)), sim_options.Fs/sim_options.Inter);
    subplot(4,1,4);
    snr(x_after_adc_int(1:length(x_after_adc_int)), sim_options.Fs/sim_options.Inter);

end
