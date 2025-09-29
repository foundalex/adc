function [data_outd, data_out] = filter_transversal(input_data_double, coeff_double, input_data, coeff, width)

    mult1 = (input_data(1) * coeff(1)); % fi(1,12,11) * fi(1,32,11)
    mult1 = fi(bitshift(mult1, -9),1,28,0);

    mult2 = (input_data(2) * coeff(2)); 
    mult2 = fi(bitshift(mult2, -9),1,28,0);

    mult3 = (input_data(3) * coeff(3)); 
    mult3 = fi(bitshift(mult3, -9),1,28,0);

    mult4 = (input_data(4) * coeff(4)); 
    mult4 = fi(bitshift(mult4, -9),1,28,0);

    mult5 = (input_data(5) * coeff(5)); 
    mult5 = fi(bitshift(mult5, -9),1,28,0);

    % data_out = fi(mult1 + mult2 + mult3 + mult4 + mult5,1,(width+11+1),0); % fi(1,86,66)
    % data_out = fi(bitshift(data_out, -width),1,12,0); % fi(1,25,0)

    data_out = (mult1 + mult2 + mult3 + mult4 + mult5); % fi(1,86,66)
    data_out = fi(bitshift(data_out,-2),1,12,0); % fi(1,28,0)

    %% double

    mult1d = input_data_double(1) * double(coeff(1)); % fi(1,82,66)
    mult2d = input_data_double(2) * double(coeff(2)); 
    mult3d = input_data_double(3) * double(coeff(3)); 
    mult4d = input_data_double(4) * double(coeff(4)); 
    mult5d = input_data_double(5) * double(coeff(5)); 

    data_outd = mult1d + mult2d + mult3d + mult4d + mult5d;



end