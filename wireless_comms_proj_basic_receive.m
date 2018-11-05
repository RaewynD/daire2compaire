clear
close all
clc

load receivedsignal.mat
load transmitsignal.mat

x_received = 10*receivedsignal;
x_transmitted = transmitsignal;

%% Grab and separate into two files

yI = real(x_transmitted);
yQ = imag(x_transmitted);

%% Define User Values
qam = 4;
L = 5;
graph = 1;
D = 1;
LL = 200; % Total number of bits
T = 1; % Symbol period in microsec
N = 21; % length of filter in symbol periods
alpha = 0.2; % alpha of sqrt raised cosine filter
fc = 5; % Carrier Frequency in MHz
fs = 200; % Sampling frequency in MHz

%% Define Pulse

p = firrcos(Ns,1/2/T,alpha,fs,'rolloff','sqrt'); 
p = p/norm(p)/sqrt(1/fs); % '1/fs' simply serves as 'delta' to approximate integral as sum

%% Define Matched Filter



%Part A - define constellation
qam_range = 1:sqrt(qam);
qam_range = d*qam_range - 0.5*d - sqrt(qam)/2;
constellation = [];
 
for xi = qam_range
    for xq = qam_range
        constellation = [constellation, xi + 1i*xq];
    end
end

% Part C - make decision
x_hat = zeros(L,1);

P_error = 0;
for x_count = 1:L
    [min_dist, index] = min(abs(constellation - x_received(x_count)));
    x_hat(x_count) = constellation(index);
    if(x_hat(x_count) ~= x_transmitted(x_count))
        P_error = P_error + 1;
    end
end

% Calculate and print error performance
Pe = P_error / L;
fprintf('Symbols wrong: %f Symbol error rate Pe = %f\n',P_error, Pe)

if graph == 1
    constellationmarkersize = 6;

    close all
    figure(1)
    LargeFigure(gcf, 0.15); % Make figure large
    clf
    zoom off
    plot(constellation,'rs','MarkerSize',constellationmarkersize,'MarkerFaceColor','r')
    set(gca,'DataAspectRatio',[1 1 1])
    grid on
    hold on
    %D = max(sqrt(Ex)*1.5, sigma_n*1.5);
    axis([-D D -D D])
    plot([-D:D/100:D],zeros(size([-D:D/100:D])),'k','LineWidth',2)
    plot(zeros(size([-D:D/100:D])),[-D:D/100:D],'k','LineWidth',2)
    xlabel('x^I, z^I')
    ylabel('x^Q, z^Q')

    title('Constellation Plot')
    plot(constellation,'rs','MarkerSize',constellationmarkersize,'MarkerFaceColor','r')
    for ii=1:L
        plot(x_received(ii),'bx')
        plot(constellation,'rs','MarkerSize',constellationmarkersize,'MarkerFaceColor','r')
        if (rem(ii,100)==0)
            pause(.00002)
        end
    end
end
