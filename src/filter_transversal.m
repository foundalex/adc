function [data_outd, data_out] = filter_transversal(input_data_double, coeff_double, input_data, coeff, width)

    mult1 = (input_data(1) * coeff(1)); % fi(1,82,66)
    mult2 = (input_data(2) * coeff(2)); 
    mult3 = (input_data(3) * coeff(3)); 
    mult4 = (input_data(4) * coeff(4)); 
    mult5 = (input_data(5) * coeff(5)); 

    % data_out = (mult1 + mult2 + mult3 + mult4 + mult5); % fi(1,86,66)


    data_out = fi(mult1 + mult2 + mult3 + mult4 + mult5,1,(width+11+1),0); % fi(1,86,66)
    data_out = fi(bitshift(data_out, -width),1,12,0);

    % if (double(data_out)) >= (2^(width+11))-1
    %     data_out = fi((2^27)-1, 1,86,0);
    % elseif (double(data_out)) <= -2^(width+11)
    %     data_out = fi(-2^27, 1,86,0);
    % end


    mult1d = input_data_double(1) * coeff_double(1); % fi(1,82,66)
    mult2d = input_data_double(2) * coeff_double(2); 
    mult3d = input_data_double(3) * coeff_double(3); 
    mult4d = input_data_double(4) * coeff_double(4); 
    mult5d = input_data_double(5) * coeff_double(5); 

    data_outd = mult1d + mult2d + mult3d + mult4d + mult5d;

    if data_outd >= 1
        data_outd = 1;
    elseif data_outd <= -1
        data_outd = -1;
    end

end