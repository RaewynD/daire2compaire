

clear
close all
clc

load receivedsignal.mat
load transmitsignal.mat

x_received = 10*receivedsignal;
x_transmitted = transmitsignal;

%user defined values
qam = 4;
d = 1;
L = 5;
graph = 1;
D = 1;

%Part A - define constellation
qam_range = 1:sqrt(qam);
qam_range = d*qam_range - 0.5*d - sqrt(qam)/2;
constellation = [];
 
for xi = qam_range
    for xq = qam_range
        constellation = [constellation, xi + j*xq];
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

    title('Scatter plot')
    pause;

    title(['Constellation Plot')
    plot(constellation,'rs','MarkerSize',constellationmarkersize,'MarkerFaceColor','r')
    for (ii=1:L)
        plot(x_received(ii),'bx')
        plot(constellation,'rs','MarkerSize',constellationmarkersize,'MarkerFaceColor','r')
        if (rem(ii,100)==0)
            pause(.00002)
        end
    end
end
