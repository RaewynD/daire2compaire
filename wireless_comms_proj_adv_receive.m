%% Emmanuel Aire-Oaihimire and Raewyn Duvall
%  Team: Daire2Compaire
%  18-758 Wireless Communications
%  Fall 2018
 
%% --Main Receive Code-- %%

clear
close all
%clc
rng('default');

% Define User Values
rake = 0; % set rake to 1 to have adding right after timing, 0 to have adding with de-spreading
real_time = 0; % set to 1 for real time, 0 for AWGN
if real_time == 1
    rot = pi/4;
    max_min = 1e5;
    freq_est_start = 100;
    freq_est_end = freq_est_start+1000;
else
    rot = 7*pi/4;
    max_min = 2e5;
    freq_est_start = 1;
    freq_est_end = freq_est_start+1000;
end
noise = 10;

load transmitsignal.mat

if real_time == 1
    load receivedsignal.mat
else
    M = 4; % M-QAM
    d = 1; % Minimum distance 
    SNR_mfb_dB = noise; % SNR_MFB in dB.  
    E_x = d^2/(6*(M-1)); % Calculate the Symbol Energy
    SNR_mfb_dB = 10^(SNR_mfb_dB/10); % Calculate the SNR
    sigma = sqrt(E_x/SNR_mfb_dB); % Calculate the STD Dev
    trans = [randn(floor(randn(1)*10),1);transmitsignal];
    receivedsignal = trans + sigma/sqrt(2)*(randn(size(trans))+j*randn(size(trans)));
end

receivedsignal = trans;

global spreading_gain spreading_mask timing_spread pilot_spread

load global_vars
% d fs Ts fc Tc T_sym F_sym symLen a p timing pilot msg ...
%    N Ns num_msg pilot_plot bits imdim msg_size ...
%    spreading_gain spreading_mask timing_spread pilot_spread_len msg_spread_len

qam = 4;
D = 1;
lenp = length(p);
L = length(msg);

% Matched filter
w = flipud(p);

x_transmitted = transmitsignal;
y_received = receivedsignal;

%% Plot

% graph = 0;
% 
% figure(8)
% LargeFigure(gcf, 0.15); % Make figure large
% clf
% 
% subplot(3,1,1)
% plot(w)
% title('Matched Filter')
% set(gca,'fontsize', 15)
% 
% s1(1) = subplot(3,1,2)
% plot(real(x_transmitted),'b')
% hold on;
% plot(imag(x_transmitted),'r')
% plot(imag(pilot_plot),'y')
% plot(real(pilot_plot),'g')
% title('Transmitted Signal')
% set(gca,'fontsize', 15)
% 
% s1(2) = subplot(3,1,3)
% plot(real(y_received),'Color',[0,0,0.7])
% hold on
% plot(imag(y_received),'Color',[0,0,0.5])
% 
% title('Received Signal')
% set(gca,'fontsize', 15)

graph = 1;

J = 2;

if graph == 1
    
    figure(J)
    LargeFigure(gcf, 0.15); % Make figure large
    clf
    subplot(3,2,1); %Filter
    plot([1:length(p)]/fs,w)
    ylabel('$p^{transmit}(t)$')
    title('Pulse Signal p(t)')
    set(gca,'fontsize', 15)
    subplot(3,2,2); %Frequency response of filter
    plot([-lenp/2+1:lenp/2]/lenp*fs,20*log10(abs(fftshift(1/sqrt(lenp)*fft(w)))))
    ylabel('$|P^{transmit}(f)|$')
    xlabel('Frequency')
    title('Frequency Response of Pulse')
    axis([-4*fc 4*fc -Inf Inf])
    set(gca,'fontsize', 15)
    subplot(3,2,3); %Transmitted Signal
    plot(real(x_transmitted),'b')
    hold on;
    plot(imag(x_transmitted),'r')
    plot(imag(pilot_plot),'y')
    plot(real(pilot_plot),'g')
    ylabel('$x^{I}(t)$, $x^{Q}(t)$')
    xlabel('Time in samples')
    legend('real','imag')
    title('Transmitted Signal')
    set(gca,'fontsize', 15)
    subplot(3,2,4); %Signal Frequency Response
    plot([0:length(transmitsignal)-1]/length(transmitsignal)-0.5, abs(fftshift(fft(transmitsignal))))
    ylabel('$|X^{base}(f)|$')
    xlabel('Frequency in 1/samples')
    title('Transimt Signal in Frequency Domain')
    set(gca,'fontsize', 15)
    subplot(3,2,5);
    plot(real(y_received),'b')
    hold on;
    plot(imag(y_received),'r')
    legend('real','imag')
    ylabel('$y^{I}(t)$,  $y^{Q}(t)$')
    xlabel('Time in samples')
    title('Received Signal')
    set(gca,'fontsize', 15)
    subplot(3,2,6)
    plot([0:length(receivedsignal)-1]/length(receivedsignal)-0.5, abs(fftshift(fft(receivedsignal))))
    ylabel('$|Y^{base}(f)|$')
    xlabel('Frequency in 1/samples')
    title('Received Signal in Frequency Domain')
    set(gca,'fontsize', 15)

    J = J+1;

end


%linkaxes(s1,'x');

%% --Greeting to the user-- %%

disp(' ')
disp('Hello. Thanks for running this program. I will try to find your message now')
disp(' ')

%% --Apply Frequency Sync-- %%

% figure(1)
% hold on;
% freq_est = y_received(freq_est_start : freq_est_end);
% freq_dom = abs(fftshift(fft(freq_est)));
% plot([0:length(freq_est)-1]/length(freq_est)-0.5, freq_dom, 'b')
% 
% title('Frequency');
% set(gca,'fontsize',15);


if graph == 1

    figure(J)
    LargeFigure(gcf, 0.15); % Make figure large
    clf
    hold on;
    grid on;
    box on;
    freq_est = y_received(freq_est_start : freq_est_end);
    freq_dom = abs(fftshift(fft(freq_est)));
    plot([0:length(freq_est)-1]/length(freq_est)-0.5, freq_dom, 'b')
    title('Frequency Response from DTFT to estimate \Delta');
    set(gca,'fontsize',15);
    
    J = J + 1;
    
end


[~, freq] = max(freq_dom);

delta_hat = (1/Ts) * freq;

%y_received = exp(-j*2*pi*delta_hat*Ts) * y_received;

%% --Apply Timing Recovery-- %%

timing_sent = timing_spread;

timing_sent = upsample(timing_sent, fs/F_sym);
timing_sent = conv(timing_sent, p);

timing_sent = reshape(timing_sent, [], 1);

[corr_time, corr_tau_time] = xcorr(timing_sent,y_received);
%[max1, offset_time1] = max(abs(corr_time1));
%[max2, offset_time2] = max(abs(corr_time1(abs(corr_time1)<max1)));
% figure()
% findpeaks(abs(corr_time));

if graph == 1
    
    figure(J)
    LargeFigure(gcf, 0.15); % Make figure large
    clf
    findpeaks(abs(corr_time));
    box on;
    grid on;
    ylabel('$y^{base}(t)$')
    title('Peak detection for Timing Recovery to create RAKE')
    set(gca,'fontsize',15);
    
    J = J+1;
    
end
    

%% --Rake-- %%

[pks, locs] = findpeaks(abs(corr_time),'SortStr','descend');
max1 = pks(1);
max1_loc = locs(1);
max2 = pks(2);
max2_loc = locs(2);
max3 = pks(3);
max3_loc = locs(3);
max4 = pks(4);
max4_loc = locs(4);

num_max = 1;
if max2 > max_min
    tau_time2 = abs(corr_tau_time(max2_loc))+1;
    y_received_timing2 = y_received(tau_time2:end);
    num_max = 2;
end
if max3 > max_min
    tau_time3 = abs(corr_tau_time(max3_loc))+1;
    y_received_timing3 = y_received(tau_time3:end);
    num_max = 3;
end
if max4 > max_min
    tau_time4 = abs(corr_tau_time(max4_loc))+1;
    y_received_timing4 = y_received(tau_time4:end);
    num_max = 3;
end

tau_time1 = abs(corr_tau_time(max1_loc))+1;
y_received_timing1 = y_received(tau_time1:end);

% tau are the actual offsets
% corr tau = offsets of correlations

if num_max == 4
    len = min([length(y_received_timing1),length(y_received_timing2),length(y_received_timing3),length(y_received_timing4)]);
    if rake == 1
        y_received_timing = y_received_timing1(1:len) + y_received_timing2(1:len) + y_received_timing3(1:len)+ y_received_timing4(1:len);
    else
        y_received_timing = y_received_timing1(1:len);
        y_received_timing2 = y_received_timing2(1:len);
        y_received_timing3 = y_received_timing3(1:len);
        y_received_timing4 = y_received_timing4(1:len);
    end
elseif num_max == 3
    len = min([length(y_received_timing1),length(y_received_timing2),length(y_received_timing3)]);
    if rake == 1
        y_received_timing = y_received_timing1(1:len) + y_received_timing2(1:len) + y_received_timing3(1:len);
    else
        y_received_timing = y_received_timing1(1:len);
        y_received_timing2 = y_received_timing2(1:len);
        y_received_timing3 = y_received_timing3(1:len);
    end
elseif num_max == 2
    len = min([length(y_received_timing1),length(y_received_timing2)]);
    if rake == 1
        y_received_timing = y_received_timing1(1:len) + y_received_timing2(1:len);
    else
        y_received_timing = y_received_timing1(1:len);
        y_received_timing2 = y_received_timing2(1:len);
    end
else
    y_received_timing = y_received_timing1;
end
    
%% Plot

% figure(10)
% LargeFigure(gcf, 0.15); % Make figure large
% clf
% subplot(3,2,1)
% stem(timing)
% title('Timing Signal')
% set(gca,'fontsize', 15)
% subplot(3,2,2)
% scatter(real(timing_sent),imag(timing_sent))
% title('Quantized Timing')
% set(gca,'fontsize', 15)
% subplot(3,2,5)
% plot(real(y_received),'b')
% hold on;
% plot(imag(y_received),'r')
% title('Received Signal')
% set(gca,'fontsize', 15)
% subplot(3,2,3)
% plot(abs(corr_time),'b')
% %hold on
% %plot(abs(corr_time2),'r')
% %plot(abs(corr_time3),'g')
% %plot(abs(corr_time4),'y')
% title('Time Correlation (Time)')
% set(gca,'fontsize', 15)
% subplot(3,2,4)
% plot(corr_tau_time,'b')
% %hold on
% %plot(corr_tau_time2,'r')
% %plot(corr_tau_time3,'g')
% %plot(corr_tau_time4,'y')
% title('Time Correlation (Time Tau)')
% set(gca,'fontsize', 15)
% subplot(3,2,6)
% plot(real(y_received_timing),'b')
% hold on;
% plot(imag(y_received_timing),'r')
% title('Y - Time Recovered')
% set(gca,'fontsize', 15)

if graph == 1

    figure(J)
    LargeFigure(gcf, 0.15); % Make figure large
    clf
    h(1) = subplot(3,2,1);
    scatter(real(timing_sent),imag(timing_sent))
    box on;
    line([0 0], ylim)
    line(xlim, [0 0])
    ylabel('$Timing^Q(t)$')
    xlabel('$Timing^I(t)$')
    title('Timing Signal - Complexspace')
    set(gca,'fontsize', 15)
    h(2) = subplot(3,2,2);
    plot(abs(corr_time))
    box on;
    ylabel('$y^{base}(t)$')
    title('Absolute Value Plot of the Time Correlation')
    set(gca,'fontsize', 15)
    h(3) = subplot(3,2,3);
    plot(corr_tau_time)
    box on;
    title('Time Correlation \tau')
    set(gca,'fontsize', 15)
    h(4) = subplot(3,2,4);
    scatter(real(y_received_timing),imag(y_received_timing))
    ylabel('$y^Q(t)$')
    xlabel('$y^I(t)$')
    box on;
    line([0 0], ylim)
    line(xlim, [0 0])
    title('$y(t)$ Recovered - Complexspace')
    set(gca,'fontsize', 15)
    h(5) = subplot(3,2,5);
    plot(real(y_received_timing),'b')
    hold on;
    plot(imag(y_received_timing),'r')
    pos = get(h,'Position');
    new = mean(cellfun(@(v)v(1),pos(1:2)));
    set(h(5),'Position',[new,pos{end}(2:end)])
    box on;
    zoom on;
    legend('real','imag')
    ylabel('Received Signal')
    xlabel('Time in samples')
    title('$y(t)$ Time Recovered Signal')
    set(gca,'fontsize', 15)

    J = J+1;

end


y_received = y_received_timing;

%% --Grab and separate into REAL and IMAGINARY-- %%

if rake == 0
    if num_max > 3
        yI4 = real(y_received_timing4);
        yQ4 = imag(y_received_timing4);
    end
    if num_max > 2
        yI3 = real(y_received_timing3);
        yQ3 = imag(y_received_timing3);
    end
    if num_max > 1
        yI2 = real(y_received_timing2);
        yQ2 = imag(y_received_timing2);
    end
end

yI = real(y_received_timing);
yQ = imag(y_received_timing);

%% Plot

% 
% x = tau_time:length(y_received);
% % x=x';
% figure(20); 
% LargeFigure(gcf, 0.15); % Make figure large
% clf
% plot(real(y_received),'b'); 
% hold on;
% plot(imag(y_received),'r');
% plot(x,yI,'g');
% plot(x,yQ,'y');
% zoom on;
% title('Overlay of Y_received and Y_correlated')
% set(gca,'fontsize', 15)

%% --Filter low pass signals with matched filter in each arm-- %%

% '1/fs' simply serves as 'delta' to approximate integral as sum

if rake == 0
    if num_max > 3
        zI4 = conv(w,yI4)*(1/fs);
        zQ4 = conv(w,yQ4)*(1/fs);
        z4_test = conv(w,y_received_timing4)*(1/fs);
    end
    if num_max > 2
        zI3 = conv(w,yI3)*(1/fs);
        zQ3 = conv(w,yQ3)*(1/fs);
        z3_test = conv(w,y_received_timing3)*(1/fs);
    end
    if num_max > 1
        zI2 = conv(w,yI2)*(1/fs);
        zQ2 = conv(w,yQ2)*(1/fs);
        z2_test = conv(w,y_received_timing2)*(1/fs);
    end
end

zI = conv(w,yI)*(1/fs);
zQ = conv(w,yQ)*(1/fs);
z_test = conv(w,y_received_timing)*(1/fs);

%% Plot

% figure(11)
% LargeFigure(gcf, 0.15); % Make figure large
% clf
% subplot(2,1,1)
% plot(zI,'b')
% hold on;
% plot(zQ,'r')
% title('Post LPF')
% set(gca,'fontsize', 15)
% subplot(2,1,2)
% scatter(zI,zQ)
% box on;
% grid on;
% hold on;
% line([0 0], ylim)
% line(xlim, [0 0])
% title('Post LPF in bitspace')
% set(gca,'fontsize', 15)

if graph == 1
    
    figure(J)
    LargeFigure(gcf, 0.15); % Make figure large
    clf
    subplot(2,1,1)
    plot(zI,'b')
    hold on;
    plot(zQ,'r')
    ylabel('Received Signal')
    xlabel('Time in samples')
    legend('real','imag')
    title('y_{base}(t) - Post Matched Filter')
    set(gca,'fontsize', 15)
    subplot(2,1,2)
    scatter(zI,zQ)
    box on;
    grid on;
    hold on;
    line([0 0], ylim)
    line(xlim, [0 0])
    ylabel('$z^Q(t)$')
    xlabel('$z^I(t)$')
    title('y_{base}(t) - Post Matched Filter in Complexspace')
    set(gca,'fontsize', 15)

    J = J+1;

end


%% --Sample filtered signal - starts falling apart here-- %%

if rake == 0
    if num_max > 3
        zIk_with_timing4 = zI4(length(w):fs*T_sym:end); 
        zQk_with_timing4 = zQ4(length(w):fs*T_sym:end);
        zk_with_timing_test4 = z4_test(length(w):fs*T_sym:end);
    end
    if num_max > 2
        zIk_with_timing3 = zI3(length(w):fs*T_sym:end); 
        zQk_with_timing3 = zQ3(length(w):fs*T_sym:end);
        zk_with_timing_test3 = z3_test(length(w):fs*T_sym:end);
    end
    if num_max > 1
        zIk_with_timing2 = zI2(length(w):fs*T_sym:end); 
        zQk_with_timing2 = zQ2(length(w):fs*T_sym:end);
        zk_with_timing_test2 = z2_test(length(w):fs*T_sym:end);
    end
end

zIk_with_timing = zI(length(w):fs*T_sym:end); 
zQk_with_timing = zQ(length(w):fs*T_sym:end);
zk_with_timing_test = z_test(length(w):fs*T_sym:end);

%% Plot

zIk_with_timing_up = upsample(zIk_with_timing,fs/F_sym);
zIk_with_timing_up_plot = [zeros(length(w)-1, 1); zIk_with_timing_up];

zQk_with_timing_up = upsample(zQk_with_timing,fs/F_sym);
zQk_with_timing_up_plot = [zeros(length(w)-1, 1); zQk_with_timing_up];

timingI_sent_plot = conv(w,real(timing_sent));
timingQ_sent_plot = conv(w,imag(timing_sent));

len_timingI = length(timingI_sent_plot) - length(w);
zI_sans_timing = zI(len_timingI : end);
zI_sans_timing_plot = [zeros(len_timingI-1, 1); zI_sans_timing];

len_timingQ = length(timingQ_sent_plot) - length(w);
zQ_sans_timing = zQ(len_timingQ : end);
zQ_sans_timing_plot = [zeros(len_timingQ-1, 1); zQ_sans_timing];

% figure(12)
% LargeFigure(gcf, 0.15); % Make figure large
% clf
% zz(1) = subplot(3,1,1);
% zoom on;
% plot(zI,'b') 
% hold on;
% stem(zIk_with_timing_up_plot,'r') 
% plot(timingI_sent_plot./10e11,'g')
% plot(zI_sans_timing_plot,'y')
% title('Post LPF signal (zI) and sampled (zIk)')
% set(gca,'fontsize', 15)
% zz(2) = subplot(3,1,2);
% zoom on;
% plot(zQ,'b') 
% hold on; 
% stem(zQk_with_timing_up_plot,'r') 
% plot(timingQ_sent_plot./10e11,'g')
% plot(zQ_sans_timing_plot,'y')
% title('Post LPF signal (zQ) and sampled (zQk)')
% set(gca,'fontsize', 15)
% subplot(3,1,3)
% scatter(zIk_with_timing,zQk_with_timing)
% box on;
% grid on;
% hold on;
% line([0 0], ylim)
% line(xlim, [0 0])
% title('Complexspace of zk')
% set(gca,'fontsize', 15)
% linkaxes(zz,'x')

if graph == 1

    figure(J)
    LargeFigure(gcf, 0.15); % Make figure large
    clf
    zz(1) = subplot(3,1,1);
    zoom on;
    plot(zI,'b') 
    hold on;
    stem(zIk_with_timing_up_plot,'r') 
    plot(timingI_sent_plot./10e11,'g')
    plot(zI_sans_timing_plot,'y')
    title('z^I and sampled z_k^I')
    legend('z^I','z_k^I - no timing','timing signal - real','z_k^I - with timing','location','northeastoutside')
    set(gca,'fontsize', 15)
    zz(2) = subplot(3,1,2);
    zoom on;
    plot(zQ,'b') 
    hold on; 
    stem(zQk_with_timing_up_plot,'r') 
    plot(timingQ_sent_plot./10e11,'g')
    plot(zQ_sans_timing_plot,'y')
    title('z^Q and sampled z_k^Q')
    legend('z^Q','z_k^Q - no timing','timing signal - imag','z_k^Q - with timing','location','northeastoutside')
    set(gca,'fontsize', 15)
    subplot(3,1,3)
    scatter(zIk_with_timing,zQk_with_timing)
    box on;
    grid on;
    hold on;
    line([0 0], ylim)
    line(xlim, [0 0])
    ylabel('$z_k^Q$')
    xlabel('$z_k^I$')
    title('Complexspace of zk')
    set(gca,'fontsize', 15)
    linkaxes(zz,'x')

    J = J+1;

end

%% --Remove Timing Preamble-- %%

if rake == 0
    if num_max > 3
        zIk_sans_timing4 = zIk_with_timing4(length(timing_spread) + 1 : end);
        zQk_sans_timing4 = zQk_with_timing4(length(timing_spread) + 1: end);
        zk_sans_timing_test4 = zk_with_timing_test4(length(timing_spread) + 1 : end);

        zk4 = zIk_sans_timing4 + j * zQk_sans_timing4;
    end
    if num_max > 2
        zIk_sans_timing3 = zIk_with_timing3(length(timing_spread) + 1 : end);
        zQk_sans_timing3 = zQk_with_timing3(length(timing_spread) + 1: end);
        zk_sans_timing_test3 = zk_with_timing_test3(length(timing_spread) + 1 : end);

        zk3 = zIk_sans_timing3 + j * zQk_sans_timing3;
    end
    if num_max > 1
        zIk_sans_timing2 = zIk_with_timing2(length(timing_spread) + 1 : end);
        zQk_sans_timing2 = zQk_with_timing2(length(timing_spread) + 1: end);
        zk_sans_timing_test2 = zk_with_timing_test2(length(timing_spread) + 1 : end);

        zk2 = zIk_sans_timing2 + j * zQk_sans_timing2;
    end
end

zIk_sans_timing = zIk_with_timing(timing_spread_len + 1 : end);
zQk_sans_timing = zQk_with_timing(timing_spread_len + 1 : end);
zk_sans_timing_test = zk_with_timing_test(timing_spread_len + 1 : end);

zk = zIk_sans_timing + j * zQk_sans_timing;

msg_mask_eq = spreading_mask(freq_spread_len + timing_spread_len + 1 : end);

%% Plot

zIk_sans_timing_up = upsample(zIk_sans_timing,fs/F_sym);
zQk_sans_timing_up = upsample(zQk_sans_timing,fs/F_sym);

% figure(13)
% LargeFigure(gcf, 0.15); % Make figure large
% clf
% zy(1) = subplot(3,1,1);
% zoom on;
% plot(zI_sans_timing,'b') 
% hold on; 
% stem(zIk_sans_timing_up,'r') 
% title('New Timing signal (zI) and sampled (zIk)')
% set(gca,'fontsize', 15)
% zy(2) = subplot(3,1,2);
% zoom on;
% plot(zQ,'b') 
% hold on; 
% stem(zQk_sans_timing_up,'r') 
% title('New Timing signal (zQ) and sampled (zQk)')
% set(gca,'fontsize', 15)
% subplot(3,1,3)
% scatter(real(zk),imag(zk))
% box on;
% grid on;
% line([0 0], ylim)
% line(xlim, [0 0])
% title('Complexspace of zk')
% set(gca,'fontsize', 15)
% linkaxes(zy,'x')

%% --Equalization and De-Spread-- %%

if rake == 0
    if num_max > 3
        zk_orig4 = zk4;

        vk_all_despread4 = [];
        vk_all4 = [];
        ho_hat_pops4 = [];
    end
    if num_max > 2
        zk_orig3 = zk3;

        vk_all_despread3 = [];
        vk_all3 = [];
        ho_hat_pops3 = [];
    end
    if num_max > 1
        zk_orig2 = zk2;

        vk_all_despread2 = [];
        vk_all2 = [];
        ho_hat_pops2 = [];
    end
end

zk_orig = zk;
msg_mask_orig = msg_mask_eq;

zk_pilot_start = 1;
zk_pilot_end = zk_pilot_start + pilot_spread_len - 1;
zk_msg_start = zk_pilot_end + 1;
zk_msg_end = zk_msg_start + msg_spread_len - 1;
zk_start = zk_msg_end + 1;

vk_all_despread = [];
vk_all = [];
ho_hat_pops = [];

if graph == 1
    
    figure(J)
    LargeFigure(gcf, 0.15); % Make figure large
    clf
    
end

for cnt = 1:num_msg
    
    msgk_pilot_mask = msg_mask_eq(zk_pilot_start : zk_pilot_end);
    msgk_mask = msg_mask_eq(zk_msg_start : zk_msg_end);
    
    % One Tap Channel
    if rake == 0
        if num_max > 3
            zk_pilot4 = zk4(zk_pilot_start : zk_pilot_end);
            zk_msg4 = zk4(zk_msg_start : zk_msg_end);
            
            ho_hat4_base = dot(msgk_pilot_mask, zk_pilot4)/norm(msgk_pilot_mask)^2;
            ho_hat4 = conj(ho_hat4_base);
            ho_hat_pops4 = [ho_hat_pops4,ho_hat4_base];
            
            %zk_msg4 = zk_msg4 * ho_hat4_base;
            
            vk_msg_spread4 = [];
            for x = 0:msg_size/2-1
                sym_start = x*spreading_gain+1;
                sym_end = (x+1)*spreading_gain;
                vk_sym_spread = zk_msg4(sym_start : sym_end) .* msgk_mask(sym_start : sym_end);
                vk_sym_despread = sum(vk_sym_spread);
                vk_msg_spread4 = [vk_msg_spread4 ; vk_sym_despread];
            end
            
            vk4 = ho_hat4 * vk_msg_spread4;
            %vk4 = vk_msg_spread4;
            
            zk4 = zk4(zk_start : end);
            disp(['The size of zk4 is: ' num2str(size(zk4))])

            vk_all_despread2 = [vk_all_despread2, vk_msg_spread2];
            vk_all4 = [vk_all4; vk4];     
        end
        if num_max > 2
            zk_pilot3 = zk3(zk_pilot_start : zk_pilot_end);
            zk_msg3 = zk3(zk_msg_start : zk_msg_end);
            
            ho_hat3_base = dot(msgk_pilot_mask, zk_pilot3)/norm(msgk_pilot_mask)^2;
            ho_hat3 = conj(ho_hat3_base);
            ho_hat_pops3 = [ho_hat_pops3,ho_hat3];
            
            %zk_msg3 = zk_msg3 * ho_hat3_base;
            
            vk_msg_spread3 = [];
            for x = 0:msg_size/2-1
                sym_start = x*spreading_gain+1;
                sym_end = (x+1)*spreading_gain;
                vk_sym_spread = zk_msg3(sym_start : sym_end) .* msgk_mask(sym_start : sym_end);
                vk_sym_despread = sum(vk_sym_spread);
                vk_msg_spread3 = [vk_msg_spread3 ; vk_sym_despread];
            end
            
            vk3 = ho_hat3 * vk_msg_spread3;
            %vk3 = vk_msg_spread3;
            
            zk3 = zk3(zk_start : end);
            disp(['The size of zk3 is: ' num2str(size(zk3))])

            vk_all_despread3 = [vk_all_despread3, vk_msg_spread3];
            vk_all3 = [vk_all3; vk3];
        end
        if num_max > 1
            zk_pilot2 = zk2(zk_pilot_start : zk_pilot_end);
            zk_msg2 = zk2(zk_msg_start : zk_msg_end);
            
            ho_hat2_base = dot(msgk_pilot_mask, zk_pilot2)/norm(msgk_pilot_mask)^2;
            ho_hat2 = conj(ho_hat2_base);
            ho_hat_pops2 = [ho_hat_pops2,ho_hat2];
            
            %zk_msg2 = zk_msg2 * ho_hat2_base;
            
            vk_msg_spread2 = [];
            for x = 0:msg_size/2-1
                sym_start = x*spreading_gain+1;
                sym_end = (x+1)*spreading_gain;
                vk_sym_spread = zk_msg2(sym_start : sym_end) .* msgk_mask(sym_start : sym_end);
                vk_sym_despread = sum(vk_sym_spread);
                vk_msg_spread2 = [vk_msg_spread2 ; vk_sym_despread];
            end
            
            vk2 = ho_hat2 * vk_msg_spread2;
            %vk2 = vk_msg_spread3;
            
            zk2 = zk2(zk_start : end);
            disp(['The size of zk2 is: ' num2str(size(zk2))])

            vk_all_despread2 = [vk_all_despread2, vk_msg_spread2];
            vk_all2 = [vk_all2; vk2];    
        end
    end

    zk_pilot = zk(zk_pilot_start : zk_pilot_end);
    zk_msg = zk(zk_msg_start : zk_msg_end);

    ho_hat_base = dot(msgk_pilot_mask, zk_pilot)/norm(msgk_pilot_mask)^2;
    ho_hat = conj(ho_hat_base);
    ho_hat_pops = [ho_hat_pops,ho_hat];
            
    %zk_msg = zk_msg * ho_hat_base;

    vk_msg_spread = [];
    for x = 0:msg_size/2-1
        sym_start = x*spreading_gain+1;
        sym_end = (x+1)*spreading_gain;
        vk_sym_spread = zk_msg(sym_start : sym_end) .* msgk_mask(sym_start : sym_end);
        vk_sym_despread = sum(vk_sym_spread);
        vk_msg_spread = [vk_msg_spread ; vk_sym_despread];
    end

    vk = ho_hat * vk_msg_spread;
    %vk = vk_msg_spread;
    
    zk = zk(zk_start : end);
    disp(['The size of zk is: ' num2str(size(zk))])

    vk_all_despread = [vk_all_despread, vk_msg_spread];
    vk_all = [vk_all; vk];
    
    msg_mask_eq = msg_mask_eq(zk_start : end);
    disp(['The size of mask is: ' num2str(size(msg_mask_eq))])

end

vk_all = vk_all / max(vk_all);
vk_all_rot = vk_all;
% vk_all_rot = vk_all * exp(j*rot);
% for x = 1:length(vk_all_rot)
%     if (((real(vk_all_rot(x)) > 0) && (imag(vk_all_rot(x)) > 0)) || ...
%             ((real(vk_all_rot(x)) < 0) && (imag(vk_all_rot(x)) < 0)))
%         vk_all_rot(x) = -vk_all_rot(x);
%     end
% end

vk_all = vk_all_rot;

    %% Plot
    
    expectedI = bits(1:2:end)*2-1;
    expectedQ = bits(2:2:end)*2-1;
    for x = 1:length(expectedI)
        expectedI(x) = expectedI(x)+rand()/4;
        expectedQ(x) = expectedQ(x)+rand()/4;
    end
    
%     figure(14)
%     cs = subplot(2,2,[1,2]);
%     scatter(real(vk_all),imag(vk_all),'b');
%     hold on;
%     scatter(real(vk_all_rot),imag(vk_all_rot),'r');
%     scatter(expectedI,expectedQ,'g');
%     box on;
%     grid on;
%     line([0 0], ylim)
%     line(xlim, [0 0])
%     title('Complexspace of vk post equalization')
%     set(gca,'fontsize', 15)

if graph == 1
    
    figure(J)
    cs = subplot(2,2,[1,2]);
    scatter(real(vk_all),imag(vk_all),'b');
    hold on;
    scatter(real(vk_all_rot),imag(vk_all_rot),'r');
    scatter(expectedI,expectedQ,'g');
    box on;
    grid on;
    line([0 0], ylim)
    line(xlim, [0 0])
    title('Complexspace of vk post equalization')
    set(gca,'fontsize', 15)
    
end


%% --Adder-- %%

if rake == 0
    if num_max == 4
        vk_hat = vk_all4 + vk_all3 + vk_all2 + vk_all;
    elseif num_max == 3
        vk_hat = vk_all3 + vk_all2 + vk_all;
    elseif num_max == 2
        vk_hat = vk_all2 + vk_all;
    else
        vk_hat = vk_all;
    end
end
        

%% --Bit Estimation-- %%

xk_hat = [];
vIk_hat = real(vk_hat);
vQk_hat = imag(vk_hat);
for x = 1:length(vk_hat)
    xk_hat = [xk_hat,vIk_hat(x),vQk_hat(x)];
end

xk_hat = sign(xk_hat);
xk_hat = (xk_hat>0);
xk_hat = reshape(xk_hat,1,length(xk_hat));

msg_hat = xk_hat(1:100);

%% Plot

% figure(14)
% subplot(2,2,[3,4])
% %stem(real(vk_all),'b')
% stem(real(vk_all_rot),'b')
% box on;
% grid on;
% hold on;
% %stem(imag(vk_all),'r')
% stem(imag(vk_all_rot)*2-1,'r')
% %stem(real(vk_all_rot),'g')
% %stem(imag(vk_all_rot),'y--')
% stem(bits(1:2:end)*2-1,'g')
% stem(bits(2:2:end)*2-1,'y--')
% legend('vIk','vQk')
% title('All of vk')

if graph == 1
    
    figure(J)
    subplot(2,2,[3,4])
    stem(real(vk_all),'b')
    box on;
    grid on;
    hold on;
    stem(imag(vk_all),'r')
    stem(bits(1:2:end)*2-1,'g')
    stem(bits(2:2:end)*2-1,'y--')
    legend('v_k^I','v_k^Q')
    title('All of v_k')
    set(gca,'fontsize', 15)
    
    J = J+2;

end

x_symb = vk_all;

pause();

%% Display Images

length(msg_hat)
msg_size
img_size = imdim(1)*imdim(2)
msg_hat_img = msg_hat(1:img_size);

length(msg_hat)
msg_size
img_size = imdim(1)*imdim(2)
msg_hat_img = msg_hat(1:img_size);

bits_err = reshape(bits,imdim);
msg_hat_err = reshape(msg_hat_img,imdim);

error_img = zeros(imdim(1),imdim(2),3);

for x = 1:imdim(1)
    for y = 1:imdim(2)
        if bits_err(x,y) ~= msg_hat_err(x,y)
            error_img(x,y,1) = 1;
            error_img(x,y,2) = 0;
            error_img(x,y,3) = 0;
        else
            if msg_hat_err(x,y) == 1
                error_img(x,y,1) = 1;
                error_img(x,y,2) = 1;
                error_img(x,y,3) = 1;
            end
        end
    end
end

%% --Additional chat with user-- %%
pause(1);
disp(' ')
disp('Thanks for waiting. Here is what I have.')
%disp(' ')
%disp(['Was this your message?: ' num2str(xk_hat)])
%disp(' ')

%% -- BER Calculation-- %%
% Compute Bit error rate (BER)
BER = mean(msg_hat_img ~= bits);
disp(' ')
disp(['Your BER is: ' num2str(BER)])
disp(' ')

%Maybe figure out how to write the text outputs onto the graph?
%Also figure out how to graph the locations of the pilots and messages in
%the large zk graph

% figure(16)
% LargeFigure(gcf, 0.15); % Make figure large
% clf
% subplot(1,3,1)
% imshow(reshape(bits,imdim))
% title('Desired Image')
% set(gca,'fontsize', 15)
% subplot(1,3,2)
% imshow(error_img)
% title('Received Image')
% xlabel(['BER is: ' num2str(BER)])
% set(gca,'fontsize', 15)

if graph == 1

    figure(J)
    LargeFigure(gcf, 0.15); % Make figure large
    clf
    subplot(1,2,1)
    imshow(reshape(bits,imdim))
    title('Desired Image')
    set(gca,'fontsize', 15)
%     subplot(1,2,2)
%     imshow(reshape(msg_hat_img,imdim))
%     title('Received Image')
%     xlabel(['BER is: ' num2str(BER)])
%     set(gca,'fontsize', 15)
    subplot(1,2,2)
    imshow(error_img)
    title('Recieved Image')
    xlabel(['BER is: ' num2str(BER)])
    set(gca,'fontsize', 15)
    pause(0.25);
    J = J+1;

end

pause();


%% --Define Constellation-- %%
qam_range = 1:sqrt(qam);
qam_range = d*qam_range - 0.5*d - sqrt(qam)/2;
constellation = [];
 
for xi = qam_range
    for xq = qam_range
        constellation = [constellation, xi + j*xq];
    end
end

% --Making some Plots-- %%
if graph == 1
    
    constellationmarkersize = 6;
    
    figure(J) %Constellation Plot Mayne
    LargeFigure(gcf, 0.15); % Make figure large
    clf
    zoom off;
    plot(constellation,'rs','MarkerSize',constellationmarkersize,'MarkerFaceColor','r')
    set(gca,'DataAspectRatio',[1 1 1])
    grid on;
    hold on;
    D = max(D, max(abs(x_symb))+1);
    axis([-D D -D D])
    plot([-D:D/100:D],zeros(size([-D:D/100:D])),'k','LineWidth',2)
    plot(zeros(size([-D:D/100:D])),[-D:D/100:D],'k','LineWidth',2)
    set(gca,'fontsize', 15)
    xlabel('$x^{I}$, $z^{I}$')
    ylabel('$x^{Q}$, $z^{Q}$')

    title('Constellation Plot')
    %plot(constellation,'rs','MarkerSize',constellationmarkersize,'MarkerFaceColor','r')
    for ii=1:length(x_symb)
        plot(x_symb(ii),'bx')
        plot(constellation,'rs','MarkerSize',constellationmarkersize,'MarkerFaceColor','r')
        if (rem(ii,100)==0)
            pause(.00002)
        end
    end

end

%% ---Helper Functions--- %%

% despread code
function [despreaded_bits] = despread_msg(msg_spread)

    global spreading_gain spreading_mask
    
    despreaded_bits = [];
    for x = 1 : length(msg_spread) / spreading_gain
        msg = msg_spread((x-1)*spreading_gain+1:x*spreading_gain) / spreading_mask;
        despreaded_bits = [despreaded_bits,msg];
    end

end

% average rake
function sum_rakes = sum_rake(rakes)
    num_rakes = size(rakes,1);
    sum_rakes = zeros(1,size(rakes,2));
    for x = 1:num_rakes
        sum_rakes = sum_rakes + rakes(x,:);
    end
end

%% To Do List:
% Split into IQ
% Apply Timing Recovery
% Apply Matched Filter
% Apply Sampling for the Zk bits
% Apply some quantizing - Rohit recommends a One-Tap Channel
% Run process to fully lay out the bits

% Phase calculation at each pilot is something to look into
% Phase offset e^jtheta or e^jtau
% Then run equalizers
% The noise at each message/pilot can be significant enough to have a large
% impact

% Rectangular pulses have some limits. It's band unlimted. The srrc is very
% effective.

% one of these should provide an AWGN channel coming back to us which is a
% good thing.