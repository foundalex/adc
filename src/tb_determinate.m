clear all;
clc;
%%
for i = 1:1000
    % dat = rand(5); % float
    dat = randi([1 100],5,5); % int
    my_func(i) = determinate(dat, dat);
    matlab_func(i) = det(dat);
    error(i) = matlab_func(i) - my_func(i); 


    [l,u] = lu(dat);
    matlab_lu(i) = prod(diag(u));

    if error(i) > 1e-11;
        disp(sprintf(['error in cycle ' num2str(i)]))
    end
end

error = error.';
my_func = my_func';
matlab_func = matlab_func';
matlab_lu = matlab_lu';