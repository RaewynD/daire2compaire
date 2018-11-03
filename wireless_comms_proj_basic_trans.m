rng('default'); % Set seed for random number generator (for repeatability of simulation)

%user defined values
qam = 4;
d = 1;

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

x2 = x1;
x_transmitted = [];
for x = 1:x_size
    if (x2(x) == 0)
        x2(x) = -1;
    end
end

for x = x_range
    if (x*2 > x_size)
        x2 = [x2,0];
    end
    x_transmitted = [x_transmitted, x2(x*2 - 1) + (j * x2(x*2))];
end

x_transmitted = 0.5 * d * x_transmitted;

transmitsignal = x_transmitted;

save('transmitsignal.mat','transmitsignal')


%% ---Helper Functions--- %%
