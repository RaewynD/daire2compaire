rng('default'); % Set seed for random number generator (for repeatability of simulation)

%user defined values
qam = 4;
d = 1;
L = 200;

%Part A - define constellation
qam_range = 1:sqrt(qam);
qam_range = d*qam_range - 0.5*d - sqrt(qam)/2;
constellation = [];
 
for xi = qam_range
    for xq = qam_range
        constellation = [constellation, xi + j*xq];
    end
end

%Part B - define points
x1 = [1,1,0,1,0,0,0,1,0,1];
x_size = size(x1);
x_size = x_size(2);
x_range = 1:ceil(x_size/2);
x_transmitted = [];

for x = x_range
    if (x*2 > x_size)
        x_tic = [x1(x*2 - 1),0];
    else
        x_tic = [x1(x*2 - 1),x1(x*2)];
    end
    x_transmitted = [x_transmitted, pnt_to_const(x_tic, d)];
end

x_transmitted

transmitsignal = x_transmitted;

save('transmitsignal.mat','transmitsignal')


%% ---Helper Functions--- %%
% translate bits to 4qam
function x_const = pnt_to_const(x_tic, d_tic)
    if x_tic(1) == 1
        if x_tic(2) == 1
            x_const = 0.5 * d_tic * (1 + 1*j);
        else
            x_const = 0.5 * d_tic * (1 - 1*j);
        end
    else
        if x_tic(2) == 1
            x_const = 0.5 * d_tic * (-1 + 1*j);
        else
            x_const = 0.5 * d_tic * (-1 - 1*j);
        end
    end
end
