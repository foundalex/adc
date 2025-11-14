function [yric, yric_int, remainder] = single_sideband(yri, y_fractional_outInt, ymi, ymi_HilbertInt, cosi, sini, num_adc, del_proc, sim_options)

        % Algorithm for working in different zones of Nyquist

        % double
        yri_cos = cosi .* yri; 
        yri_sin = sini .* ymi; 

        if (mod(sim_options.Z,2) == 0)
            yric = yri_cos + yri_sin;
        else
            yric = yri_cos - yri_sin;
        end

        %% integer
        if sim_options.Z == 2 || sim_options.Z == 3
            if (num_adc == 1) % ADC2
                yric_int = ymi_HilbertInt;
                
                remainder = sim_options.remaind_fir_hilbert(1);
                yhil_imag_int_double = double(ymi_HilbertInt)*2^-remainder;

                relative_error1 = yric./yhil_imag_int_double;
                figure(7);
                subplot(2,1,1)
                plot([yric, yhil_imag_int_double]);
                subplot(2,1,2)
                plot(relative_error1)

            elseif (num_adc == 2) % ADC3
                yric_int = -y_fractional_outInt;

                remainder = sim_options.remaind_fir_fractional(2);
                y_fractional_outInt_double = double(y_fractional_outInt)*2^-(remainder);
                relative_error2 = yric./-y_fractional_outInt_double;
                figure(8);
                subplot(2,1,1)
                plot([yric, -y_fractional_outInt_double]);
                subplot(2,1,2)
                plot(relative_error2)

            elseif (num_adc == 3) % ADC4
                yric_int = -ymi_HilbertInt;

                remainder = sim_options.remaind_fir_hilbert(3);
                yhil_imag_int_double = double(-ymi_HilbertInt)*2^-remainder;
                relative_error3 = yric./yhil_imag_int_double;
                figure(9);
                subplot(2,1,1)
                plot([yric, yhil_imag_int_double]);
                subplot(2,1,2)
                plot(relative_error3)
            end
        elseif sim_options.Z == 4
            if (num_adc == 1)
                yric_int = -y_fractional_outInt;
                remainder = sim_options.remaind_fir_fractional(1);
            elseif (num_adc == 2)
                yric_int = y_fractional_outInt;
                remainder = sim_options.remaind_fir_fractional(2);
            elseif (num_adc == 3)
                yric_int = -y_fractional_outInt;
                remainder = sim_options.remaind_fir_fractional(3);
            end
        else 
            yric_int = y_fractional_outInt;
            remainder = sim_options.remaind_fir_fractional(num_adc);
        end
        
end