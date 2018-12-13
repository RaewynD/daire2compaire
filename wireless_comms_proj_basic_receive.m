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
srrc = 1; % "1" for SRRC to be on. "0" for rectangular pulse.
real_time = 0; % "1" for Tx/Rx to be on. "0" for AWGN. (must be opposite AWGN)
AWGN = 1; % "1" for AWGN to be on. "0" for AWGN off. (must be opposite real_time)
graph = 1; % "1" for graphs to be on. "0" for them to be off.

% Loads files from the transmitter file
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
%d fs Ts fc Tc T_sym F_sym symLen a p timing pilot msg N Ns num_msg pilot_plot bits imdim msg_size
qam = 4;
D = 1;
lenp = length(p);
L = length(msg);

% Matched filter
w = flipud(p);

y_received = receivedsignal;
x_transmitted = transmitsignal;

J = 2;

if graph == 1
    
    figure(J)
    LargeFigure(gcf, 0.15); % Make figure large
    clf
    subplot(3,2,1); %Filter
    plot([1:length(p)]/fs,p)
    ylabel('$p^{transmit}(t)$')
    title('Pulse Signal p(t)')
    set(gca,'fontsize', 15)
    subplot(3,2,2); %Frequency response of filter
    plot([-lenp/2+1:lenp/2]/lenp*fs,20*log10(abs(fftshift(1/sqrt(lenp)*fft(p)))))
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

%% --Greeting to the user-- %%

disp(' ')
disp('Hello. Thanks for running this program. I will try to find your message now')
disp(' ')
disp(['I am searching for Pilot signal: ' num2str(pilot)])
disp(' ')

if J == 1

    figure(J)
    LargeFigure(gcf, 0.15); % Make figure large
    clf
    stem(pilot)
    title('Pilot Signal')
    set(gca,'fontsize', 15)

    J = J+1;

end

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
%% --Grab and separate into REAL and IMAGINARY-- %%

yI = real(y_received_timing);
yQ = imag(y_received_timing);

x = tau_time:length(y_received);
x=x';

if graph == 1

    figure(J); 
    LargeFigure(gcf, 0.15); % Make figure large
    clf
    plot(real(y_received),'b'); 
    hold on;
    plot(imag(y_received),'r');
    plot(x,yI,'g');
    plot(x,yQ,'y');
    zoom on;
    ylabel('Received Signal')
    xlabel('Time in samples')
    legend('Orig - real','Orig - imag','Sync - real','Sync - imag')
    title('Overlay of y_{received}(t) and y_{synchronized}(t) signal')
    set(gca,'fontsize', 15)

    J = J+1;

end


%% --Filter low pass signals with matched filter in each arm-- %%

% '1/fs' simply serves as 'delta' to approximate integral as sum
zI_matched = conv(w,yI)*(1/fs);
zQ_matched = conv(w,yQ)*(1/fs);

zI = zI_matched;
zQ = zQ_matched;

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
    title('Post LPF signal z^I and sampled z_k^I')
    legend('z^I','z_k^I - no timing','timing signal - real','z_k^I - with timing','location','northeastoutside')
    set(gca,'fontsize', 15)
    zz(2) = subplot(3,1,2);
    zoom on;
    plot(zQ,'b') 
    hold on; 
    stem(zQk_with_timing_up_plot,'r') 
    plot(timingQ_sent_plot./10e11,'g')
    plot(zQ_sans_timing_plot,'y')
    title('Post LPF signal z^Q and sampled z_k^Q')
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

zIk_sans_timing = zIk_with_timing(length(timing)/2 + 1 : end);
zQk_sans_timing = zQk_with_timing(length(timing)/2 + 1: end);

zk = zIk_sans_timing + j * zQk_sans_timing;

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

%% --Equalization-- %%

zk_orig = zk;

pilot = 2*pilot - 1;
complex_pilot = pilot(1:2:end) + j*pilot(2:2:end);

zk_pilot_end = length(complex_pilot);
zk_msg_start = zk_pilot_end + 1;
zk_msg_end = zk_msg_start + msg_size/2 - 1; %msg_size is in bit space
zk_start = zk_msg_end + 1;

msg_hat = [];
vk_all = [];
ho_hat_pops = [];

if graph == 1
    
    figure(J)
    LargeFigure(gcf, 0.15); % Make figure large
    clf
    
end

for cnt = 1:num_msg
    zk_pilot = zk(1 : zk_pilot_end);
    zk_msg = zk(zk_msg_start : zk_msg_end);
    zk = zk(zk_start : end);
    disp(['The size of zk is: ' num2str(size(zk))]);
    
    ho_hat = dot(complex_pilot,zk_pilot)/norm(complex_pilot)^2;
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
    
    if graph == 1
        
        figure(J)
        cs = subplot(2,2,[1,2]);
        scatter(real(zk_msg),imag(zk_msg),'r');
        box on;
        grid on;
        %hold on;
        line([0 0], ylim)
        line(xlim, [0 0])
        ylabel('$Q$')
        xlabel('$I$')
        title('Complexspace of v_k post equalization')
        legend('v_k','z_k pilot','v_k pilot','location','northeastoutside')
        set(gca,'fontsize', 15)
    
    end
    
end

x_symb = vk_all;

if graph == 1

    figure(J+1)
    LargeFigure(gcf, 0.15); % Make figure large
    clf
    scatter(real(ho_hat_pops),imag(ho_hat_pops))
    box on;
    grid on;
    line([0 0], ylim)
    line(xlim, [0 0])
    title('Equalization values of h_o')
    set(gca,'fontsize',15)
    
end

if graph == 1
    
    figure(J)
    subplot(2,2,[3,4])
    stem(real(vk_all),'b')
    box on;
    grid on;
    hold on;
    stem(imag(vk_all),'r')
    legend('v_k^I','v_k^Q')
    title('All of v_k')
    set(gca,'fontsize', 15)
    
    J = J+2;

end

%% Display Images

length(msg_hat);
msg_size;
img_size = imdim(1)*imdim(2);
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

pause();

%% -- BER Calculation-- %%
% Compute Bit error rate (BER)
BER = mean(msg_hat_img ~= bits);
disp(' ')
disp(['Your final BER is: ' num2str(BER)])
disp(' ')

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
