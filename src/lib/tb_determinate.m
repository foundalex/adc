clear all;
clc;
%%
for i = 1:1000
    % dat = rand(5); % float
    dat = randi([-2048 2047],5,5); % int

    aa = int16(dat);

    [DetM_5x5(i), DetM_5x5_int(i)] = determinate(double(dat)*2^-11, aa);

    % er_det_procent(i) = (double(DetM_5x5_int(i))*2^-49 - DetM_5x5(i)) / DetM_5x5(i);

    er_det_procent_lu(i) = DetM_5x5(i)/det(double(dat)*2^-11);
    er_det_procent(i) = ((double(DetM_5x5_int(i))*2^-49) / DetM_5x5(i));

    % [l,u] = lu(dat);
    % matlab_lu(i) = prod(diag(u));
    % 
    % if error(i) > 1e-11;
    %     disp(sprintf(['error in cycle ' num2str(i)]))
    % end
end

er_det_procent = er_det_procent.';
% my_func = my_func';
% matlab_func = matlab_func';
% matlab_lu = matlab_lu';

figure(7)
subplot(2,1,1)
plot(er_det_procent)
subplot(2,1,2)
plot(er_det_procent_lu)