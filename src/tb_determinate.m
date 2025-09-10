clear all;
clc;
%%
for i = 1:1000
    dat = randi([1 10],5,5);
    my_func = determinate(dat, 0);
    matlab_func = det(dat);
    error = abs(matlab_func - my_func);
    if error > 0.0001
        disp(sprintf(['error in cycle ' num2str(i)]))
    end
end