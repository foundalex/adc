function [yric, yric_int] = single_sideband(yri, y_fractional_outInt, ymi, ymi_HilbertInt, cosi, sini, num_adc, del_proc, sim_options)


        yhil_imag = [ymi(del_proc+1:end); zeros(del_proc,1)];
        yhil_imag_int = [ymi_HilbertInt(del_proc+1:end); zeros(del_proc,1)];
        %% Algorithm for working in different zones of Nyquist
        %%
        % double
        yri_cos = cosi .* yri; 
        % yri_sin = sini .* ymi(del_proc+1:end); 
        yri_sin = sini .* yhil_imag; 

        if (mod(sim_options.Z,2) == 0)
            yric = yri_cos + yri_sin;
        else
            yric = yri_cos - yri_sin;
        end

        yric(del_proc+1:end);


        %% integer

        yhil_imag_int_double = double(yhil_imag_int)*2^-sim_options.divide_remainder;
        y_fractional_outInt_double = double(y_fractional_outInt)*2^-(sim_options.fractional_coeff_width-1);

        if sim_options.Z == 2 || sim_options.Z == 3
            if (num_adc == 1) % ADC2
                yric_int = yhil_imag_int;

                relative_error1 = yric(del_proc+1:end)./yhil_imag_int_double(del_proc+1:end);
                figure(7);
                subplot(2,1,1)
                plot([yric(del_proc+1:end), yhil_imag_int_double(del_proc+1:end)]);
                subplot(2,1,2)
                plot(relative_error1)

            elseif (num_adc == 2) % ADC3
                yric_int = -y_fractional_outInt;

                relative_error2 = yric./-y_fractional_outInt_double;
                figure(8);
                subplot(2,1,1)
                plot([yric(del_proc+1:end), -y_fractional_outInt_double(del_proc+1:end)]);
                subplot(2,1,2)
                plot(relative_error2)

            elseif (num_adc == 3) % ADC4
                yric_int = -yhil_imag_int;

                relative_error3 = yric./-yhil_imag_int_double;
                figure(9);
                subplot(2,1,1)
                plot([yric(del_proc+1:end), -yhil_imag_int_double(del_proc+1:end)]);
                subplot(2,1,2)
                plot(relative_error3)
            end
        elseif sim_options.Z == 4
            if (num_adc == 1)
                yric_int = -y_fractional_outInt;
            elseif (num_adc == 2)
                yric_int = y_fractional_outInt;
            elseif (num_adc == 3)
                yric_int = -y_fractional_outInt;
            end
        else 
            yric_int = y_fractional_outInt;
        end

        yric_int(del_proc+1:end);

end