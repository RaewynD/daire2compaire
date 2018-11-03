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

%Part B - define L random points
x1 = rand(L,1);
x_transmitted = ceil(x1 * qam);

for x = 1:L
    x_transmitted(x) = constellation(x_transmitted(x));
end

transmitsignal = x_transmitted;

save('transmitsignal.mat','transmitsignal')
