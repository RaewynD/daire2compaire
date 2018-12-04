%% Emmanuel Aire-Oaihimire and Raewyn Duvall
%  Team: Daire2Compaire
%  18-758 Wireless Communications
%  Fall 2018
 
%% --Main Receive Code-- %%

clear
%close all
%clc
rng('default');

% Define User Values
srrc = 1;
real_time = 0;
AWGN = 1;

%load transmitsignal.mat

if srrc == 1
    %load receivedsignal_SRRC
    load transmitsignal_SRRC
elseif srrc == 0
    %load receivedsignal_RECT
    load transmitsignal_RECT
else
end

if real_time == 1
    load receivedsignal.mat
end

if AWGN == 1
    M = 4; % M-QAM
    d = 1; % Minimum distance 
    SNR_mfb_dB = 10; % SNR_MFB in dB.  
    E_x = d^2/(6*(M-1)); % Calculate the Symbol Energy
    SNR_mfb_dB = 10^(SNR_mfb_dB/10); % Calculate the SNR
    sigma = sqrt(E_x/SNR_mfb_dB); % Calculate the STD Dev
    receivedsignal = transmitsignal + sigma/sqrt(2)*(randn(size(transmitsignal))+j*randn(size(transmitsignal)));
end

load global_vars
%d fs Ts fc Tc T_sym F_sym symLen a p timing pilot msg
qam = 4;
D = 1;
lenp = length(p);
L = length(msg);


% Matched filter
w = flipud(p);

y_received = receivedsignal;
x_transmitted = transmitsignal;

graph = 0;

figure(8)
LargeFigure(gcf, 0.15); % Make figure large
clf
subplot(3,1,1)
plot(w)
title('Matched Filter')
set(gca,'fontsize', 15)
subplot(3,1,2)
plot(real(x_transmitted),'b')
hold on;
plot(imag(x_transmitted),'r')
plot(imag(pilot_plot),'y')
plot(real(pilot_plot),'g')
title('Transmitted Signal')
set(gca,'fontsize', 15)
subplot(3,1,3)
plot(real(y_received),'b')
hold on;
plot(imag(y_received),'r')
title('Received Signal')
set(gca,'fontsize', 15)

%% --Greeting to the user-- %%

disp(' ')
disp('Hello. Thanks for running this program. I will try to find your message now')
disp(' ')
disp(['I am searching for Pilot signal: ' num2str(pilot)])
disp(' ')

figure(9)
LargeFigure(gcf, 0.15); % Make figure large
clf
stem(pilot)
title('Pilot Signal')
set(gca,'fontsize', 15)

%% --Apply Timing Recovery-- %%

timing_sent = 2*timing - 1;
timing_sent = timing_sent';
timing_I = timing_sent(1:2:end)*(0.5*d);
timing_Q = timing_sent(2:2:end)*(0.5*d);

timing_I = upsample(timing_I, fs/F_sym);
timing_Q = upsample(timing_Q, fs/F_sym);
timing_I = conv(timing_I, p);
timing_Q = conv(timing_Q, p);

timing_sent = timing_I + j*timing_Q;
timing_sent = reshape(timing_sent, [], 1);

[corr_time, corr_tau_time] = xcorr(timing_sent, y_received);
[~, offset_time] = max(abs(corr_time));
tau_time = abs(corr_tau_time(offset_time))+1;
% tau are the actual offsets
% corr tau = offsets of correlations

y_received_timing = y_received(tau_time:end);

figure(10)
LargeFigure(gcf, 0.15); % Make figure large
clf
subplot(3,2,1)
stem(timing)
title('Timing Signal')
set(gca,'fontsize', 15)
subplot(3,2,2)
scatter(real(timing_sent),imag(timing_sent))
title('Quantized Timing')
set(gca,'fontsize', 15)
subplot(3,2,5)
plot(real(y_received),'b')
hold on;
plot(imag(y_received),'r')
title('Received Signal')
set(gca,'fontsize', 15)
subplot(3,2,3)
plot(abs(corr_time))
title('Time Correlation (Time)')
set(gca,'fontsize', 15)
subplot(3,2,4)
plot(corr_tau_time)
title('Time Correlation (Time Tau)')
set(gca,'fontsize', 15)
subplot(3,2,6)
plot(real(y_received_timing),'b')
hold on;
plot(imag(y_received_timing),'r')
title('Y - Time Recovered')
set(gca,'fontsize', 15)

%% --Grab and separate into REAL and IMAGINARY-- %%

yI = real(y_received_timing);
yQ = imag(y_received_timing);

x = tau_time:length(y_received);x=x';
figure(20); 
LargeFigure(gcf, 0.15); % Make figure large
clf
plot(real(y_received),'b'); 
hold on;
plot(imag(y_received),'r');
plot(x,yI,'g');
plot(x,yQ,'y');
zoom on;
title('Overlay of Y_received and Y_correlated')
set(gca,'fontsize', 15)

%% --Filter low pass signals with matched filter in each arm-- %%

% '1/fs' simply serves as 'delta' to approximate integral as sum
zI = conv(w,yI)*(1/fs);
zQ = conv(w,yQ)*(1/fs);

figure(11)
LargeFigure(gcf, 0.15); % Make figure large
clf
subplot(2,1,1)
plot(zI,'b')
hold on;
plot(zQ,'r')
title('Post LPF')
set(gca,'fontsize', 15)
subplot(2,1,2)
scatter(zI,zQ)
box on;
grid on;
hold on;
line([0 0], ylim)
line(xlim, [0 0])
title('Post LPF in bitspace')
set(gca,'fontsize', 15)

%% --Sample filtered signal - starts falling apart here-- %%

if srrc == 1
    zIk_with_timing = zI(length(w):fs*T_sym:end); 
    zQk_with_timing = zQ(length(w):fs*T_sym:end);
else
    zIk_with_timing = zI(1:fs*T_sym:end); 
    zQk_with_timing = zQ(1:fs*T_sym:end);
end

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

figure(12)
LargeFigure(gcf, 0.15); % Make figure large
clf
zz(1) = subplot(3,1,1);
zoom on;
plot(zI,'b') 
hold on;
stem(zIk_with_timing_up_plot,'r') 
plot(timingI_sent_plot./10e11,'g')
plot(zI_sans_timing_plot,'y')
title('Post LPF signal (zI) and sampled (zIk)')
set(gca,'fontsize', 15)
zz(2) = subplot(3,1,2);
zoom on;
plot(zQ,'b') 
hold on; 
stem(zQk_with_timing_up_plot,'r') 
plot(timingQ_sent_plot./10e11,'g')
plot(zQ_sans_timing_plot,'y')
title('Post LPF signal (zQ) and sampled (zQk)')
set(gca,'fontsize', 15)
subplot(3,1,3)
scatter(zIk,zQk)
box on;
grid on;
hold on;
line([0 0], ylim)
line(xlim, [0 0])
title('Complexspace of zk')
set(gca,'fontsize', 15)
linkaxes(zz,'x')

zk_with_timing = (zIk_with_timing + j*zQk_with_timing);

zk = zk_with_timing(length(timing)+1 : end);
zIk = real(zk);
zQk = imag(zk);

figure(13)
LargeFigure(gcf, 0.15); % Make figure large
clf
zy(1) = subplot(3,1,1);
zoom on;
plot(zI,'b') 
hold on; 
stem(upsample(real(zk_with_timing),fs/F_sym),'r') 
title('New Timing signal (zI) and sampled (zIk)')
set(gca,'fontsize', 15)
zy(2) = subplot(3,1,2);
zoom on;
plot(zQ,'b') 
hold on; 
stem(upsample(imag(zk),fs/F_sym),'r') 
title('New Timing signal (zQ) and sampled (zQk)')
set(gca,'fontsize', 15)
subplot(3,1,3)
scatter(real(zk),imag(zk))
box on;
grid on;
line([0 0], ylim)
line(xlim, [0 0])
title('Complexspace of zk')
set(gca,'fontsize', 15)
linkaxes(zy,'x')

%% --Equalization-- %%

zk_orig = zk;

pilot = 2*pilot - 1;
complex_pilot = pilot(1:2:end) + j*pilot(2:2:end);
zk_pilot_end = length(complex_pilot);
zk_msg_start = zk_pilot_end + 1;
zk_msg_end = zk_msg_start + msg_size/2 - 1; %msg_size is in bit space
zk_start = zk_msg_end + 1;

msg_hat = [];
ho_hat_pops = [];
vk_all = [];


for cnt = 1:num_msg
    zk_pilot = zk(1 : zk_pilot_end);
    zk_msg = zk(zk_msg_start : zk_msg_end);
    zk = zk(zk_start : end);
    disp(['The size of zk is: ' num2str(size(zk))])
    
    %ho_hat = (w'*kz(1:length(complex_pilot)))/(w'*complex_pilot)
    %ho_hat = (complex_pilot' .* kz(1:length(complex_pilot)))/(norm(complex_pilot)^2);
    ho_hat = dot(conj(complex_pilot),zk_pilot)/norm(complex_pilot)^2;
    ho_hat_pops = [ho_hat_pops,ho_hat];

    %% Detect bits - One Tap Channel

    vk_pilot = zk_pilot / ho_hat;
    vk = zk_msg / ho_hat;
    
    vIk = real(vk);
    vQk = imag(vk);
    
    vk_bits = [];
    
    for x = 1:length(vIk)
        vk_bits = [vk_bits,vIk(x),vQk(x)];
    end
    
    xk_hat = sign(vk_bits);
    xk_hat = (xk_hat>0);
    xk_hat = reshape(xk_hat,1,length(xk_hat));
    
    msg_hat = [msg_hat,xk_hat];
    vk_all = [vk_all;vk];
    
    f = figure(14)
    if cnt == 1
        LargeFigure(gcf, 0.15); % Make figure large
        clf
    end
    %cs = subplot(2,2,[1,2]);
    cla(f);
    scatter(real(vk),imag(vk));
    hold on;
    scatter(real(zk_pilot), imag(zk_pilot));
    scatter(real(vk_pilot), imag(vk_pilot));
    box on;
    grid on;
    hold on;
    line([0 0], ylim)
    line(xlim, [0 0])
    title('Complexspace of vk post equalization')
    set(gca,'fontsize', 15)
    %subplot(2,2,[3,4])
    %stem(real(vk_all),'b')
    %box on;
    %grid on;
    %hold on;
    %title('vIk')
    %subplot(3,1,2)
    %stem(imag(vk_all),'r')
    %legend('vIk','vQk')
    %title('All of vk')
    %title('vQk')
    %pause(0.5);
    %pause();
    
    %figure(15)
    %if cnt == 1
    %    LargeFigure(gcf, 0.15); % Make figure large
    %    clf
    %end
    %subplot(2,1,1)
    %stem(real(complex_pilot),'b')
    %hold on;
    %stem(real(zk_pilot/ho_hat),'r')
    %legend('Complex Pilot','Sampled pilot')
    %box on;
    %grid on;
    %title('Real Pilot Competition')
    %set(gca,'fontsize', 15)
    %hold off;
    %subplot(2,1,2)
    %stem(imag(complex_pilot),'b')
    %hold on;
    %stem(imag(zk_pilot/ho_hat),'r')
    %legend('Complex Pilot','Sampled pilot')
    %box on;
    %grid on;
    %title('Imaginary Pilot Competition')
    %set(gca,'fontsize', 15)
    %hold off;
    %pause(0.125);
    %pause();

end

pause();

length(msg_hat)
msg_size
img_size = imdim(1)*imdim(2)
msg_hat_img = msg_hat(1:img_size);

figure(16)
LargeFigure(gcf, 0.15); % Make figure large
clf
subplot(1,2,1)
imshow(reshape(bits,imdim))
title('Desired Image')
set(gca,'fontsize', 15)
subplot(1,2,2)
imshow(reshape(msg_hat_img,imdim))
title('What is good mayne?')
set(gca,'fontsize', 15)

%figure(15)

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

%% --Define Constellation-- %%
qam_range = 1:sqrt(qam);
qam_range = d*qam_range - 0.5*d - sqrt(qam)/2;
constellation = [];
 
for xi = qam_range
    for xq = qam_range
        constellation = [constellation, xi + 1i*xq];
    end
end

%% --Making some Plots-- %%
if graph == 1
    constellationmarkersize = 6;
    ax = [];
    cats = 0;
    lono = 1:length(zk_sign);

    %close all
    figure(1) %Constellation Plot Mayne
    LargeFigure(gcf, 0.15); % Make figure large
    clf
    zoom off
    plot(constellation,'rs','MarkerSize',constellationmarkersize,'MarkerFaceColor','r')
    set(gca,'DataAspectRatio',[1 1 1])
    grid on
    hold on
    D = max(D, max(abs(v_k))+1);
    axis([-D D -D D])
    plot([-D:D/100:D],zeros(size([-D:D/100:D])),'k','LineWidth',2)
    plot(zeros(size([-D:D/100:D])),[-D:D/100:D],'k','LineWidth',2)
    set(gca,'fontsize', 15)
    xlabel('$x^{I}$, $z^{I}$')
    ylabel('$x^{Q}$, $z^{Q}$')

    title('Constellation Plot')
    plot(constellation,'rs','MarkerSize',constellationmarkersize,'MarkerFaceColor','r')
    for ii=1:L/2
        plot(v_k(ii),'bx')
        plot(constellation,'rs','MarkerSize',constellationmarkersize,'MarkerFaceColor','r')
        if (rem(ii,100)==0)
            pause(.00002)
        end
    end


    % Display signals in time and frequency domain
    figure(2)
    LargeFigure(gcf, 0.15); % Make figure large
    clf
    subplot(3,2,1);
    plot([1:lenp]/fs,p)
    ylabel('$p^{transmit}(t)$')
    title('Pulse Signal p(t)')
    set(gca,'fontsize', 15)
    subplot(3,2,3);
    plot(real(x_transmitted),'b')
    hold on
    plot(imag(x_transmitted),'r')
    %plot(imag(pilot_plot),'y')
    %plot(real(pilot_plot),'g')
    legend('real','imag')
    ylabel('$x^{I}(t)$, $x^{Q}(t)$')
    xlabel('Time in samples')
    title('Transmit Signal')
    set(gca,'fontsize', 15)
    subplot(3,2,2);
    plot([-lenp/2+1:lenp/2]/lenp*fs,20*log10(abs(fftshift(1/sqrt(lenp)*fft(p)))))
    ylabel('$|P^{transmit}(f)|$')
    title('Pulse Signal in Frequency Domain')
    %axis([-4*fc 4*fc -40 40])
    set(gca,'fontsize', 15)
    subplot(3,2,4);
    plot([0:length(transmitsignal)-1]/length(transmitsignal)-0.5, abs(fftshift(fft(transmitsignal))))
    ylabel('$|X^{base}(f)|$')
    xlabel('Frequency in 1/samples')
    title('Transimt Signal in Frequency Domain')
    set(gca,'fontsize', 15)
    subplot(3,2,5)
    plot(real(y_received),'b')
    hold on
    plot(imag(y_received),'r')
    zoom xon
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
    
    figure(3) %Timing recovery
    LargeFigure(gcf, 0.15); % Make figure large
    clf
    subplot(2,2,1);
    plot(corr_time)
    ylabel('$y^{base}(t)$')
    title('Time Correlation')
    set(gca,'fontsize', 15)
    subplot(2,2,2);
    plot(abs(corr_time))
    ylabel('$y^{base}(t)$')
    title('Absolute Value of the Time Correlation')
    set(gca,'fontsize', 15)
    subplot(2,2,3);
    %plot(y_received_timing)
    scatter(real(y_received_timing),imag(y_received_timing))
    ylabel('$y^Q(t)$')
    xlabel('$y^I(t)$')
    title('$y(t)$ Time Recovered')
    set(gca,'fontsize', 15)
    subplot(2,2,4);
    plot(real(y_received_timing),'b')
    hold on;
    plot(imag(y_received_timing),'r')
    zoom xon
    legend('real','imag')
    ylabel('Received Signal')
    xlabel('Time in samples')
    title('$y(t)$ Time Recovered Signal')
    set(gca,'fontsize', 15)
    
    figure(4) %Frame recovery
    LargeFigure(gcf, 0.15); % Make figure large
    clf
    subplot(3,2,1);
    plot(corr_frame)
    ylabel('$z_k$')
    title('Frame Correlation Start')
    set(gca,'fontsize', 15)
    subplot(3,2,2);
    plot(abs(corr_frame))
    ylabel('$z_k$')
    title('Absolute Value of the Frame Correlation Start')
    set(gca,'fontsize', 15)
    subplot(3,2,3);
    plot(corr_frame_end)
    ylabel('$z_k$')
    title('Frame Correlation End')
    set(gca,'fontsize', 15)
    subplot(3,2,4);
    plot(abs(corr_frame_end))
    ylabel('$z_k$')
    title('Absolute Value of the Frame Correlation End')
    set(gca,'fontsize', 15)
    subplot(3,2,5);
    stem(zk_bits)
    ylabel('$z_k$')
    title('Recovered $z_k$ Bits')
    set(gca,'fontsize', 15)
    subplot(3,2,6);
    stem(zk_msg)
    ylabel('$z_k$')
    title('Frame Correlation End - Found message')
    set(gca,'fontsize', 15)
    
    figure(5) % Bit Stems Plots
    LargeFigure(gcf, 0.15); % Make figure large
    clf
    ax(1) = subplot(2,2,1);
    %stem([1:x],bitI_hat,'b')
    %hold on
    stem(1:length(pilot),pilot,'b')
    %for zee = 1:length(zIk_frame)
    %    cat = find(zIk==zIk_frame(zee));
    %    stem(lono(cat),zI_bits(cat),'b')
    %end
    ylabel('Bits') %'$x^I_k,   z^I_{k}$'
    xlabel('discrete time')
    title('Pilot Signal')
    ylim([-2 2]);
    set(gca,'fontsize', 15)
    ax(2) = subplot(2,2,2);
    %stem([1:x],bitI_hat,'b')
    %hold on
    stem(lono,zk_sign,'r')
    %for zee = 1:length(zIk_frame)
    %    cat = find(zIk==zIk_frame(zee));
    %    stem(lono(cat),zI_bits(cat),'b')
    %end
    ylabel('$z_{k}$') %'$x^I_k,   z^I_{k}$'
    xlabel('discrete time  $k$  (sampled at $t=kT$)')
    title('Sampler Output $z_{k}$')
    ylim([-2 2]);
    set(gca,'fontsize', 15)
    %ax(2) = subplot(2,2,2);
    %stem([1:x],bitQ_hat,'b')
    %hold on
    %stem(lono,zQ_bits,'b')
    %for zee = 1:length(zQk_frame)
    %    cat = find(zQk==zQk_frame(zee));
    %    stem(lono(cat),zQ_bits(cat),'r')
    %end
    %ylabel('$z^Q_{k}$') %'$x^Q_k,   z^Q_{k}$'
    %xlabel('discrete time  $k$  (sampled at $t=kT$)')
    %title('Sampler Output $z^Q_{k}$')
    %ylim([-2 2]);
    %set(gca,'fontsize', 15)
    subplot(2,2,3);
    stem(1:length(xk_hat),xk_hat,'b')
    hold on
    ylabel('$v_{k}$') 
    xlabel('discrete time  $k$  (sampled at $t=kT$)')
    title('Equalizer $v_{k}$ output samples')
    set(gca,'fontsize', 15)
    ylim([-2 2]);
    %subplot(2,2,4);
    %stem([1:x],bitQ_hat,'b')
    %hold on
    %stem([1:L],bits_hat','k')
    %ylabel('$z_{k}$') %'$x^Q_k,   z^Q_{k}$'
    %xlabel('discrete time  $k$  (sampled at $t=kT$)')
    %title('Decoded Bits $z_{k}$')
    %ylim([-2 2]);
    %set(gca,'fontsize', 15)
    %linkaxes(ax,'x')
    %zoom xon

    %figure(6)
    %scatter(vIk, vQk);
    %grid on;
    
    
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