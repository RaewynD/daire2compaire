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
srrc = 2;
real_time = 0;
AWGN = 1;

load transmitsignal.mat

if srrc == 1
    load receivedsignal_SRRC
    load transmitsignal_SRRC
elseif srrc == 0
    load receivedsignal_RECT
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

graph = 1;

%% Apply Timing Recovery

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

%% Grab and separate into REAL and IMAGINARY

yI = real(y_received_timing);
yQ = imag(y_received_timing);

%% Filter low pass signals with matched filter in each arm

% '1/fs' simply serves as 'delta' to approximate integral as sum
zI = conv(w,yI)*(1/fs);
zQ = conv(w,yQ)*(1/fs);

%% Sample filtered signal - starts falling apart here

k = 1;
%zIk = [];
%zQk = [];
for s = length(w)+1:T_sym/Ts:length(zI)
    zk(k*2-1) = zI(s); %odds
    zk(k*2) = zQ(s); %evens
%    zIk(k) = zI(s);
%    zQk(k) = zQ(s);
    k = k+1;
end

%xk_should = zk(41:50);

%xkI_should = xk_should(1:2:end);
%xkQ_should = xk_should(2:2:end);

%% Frame Recovery
% remove when not doing phase recovery, too

zk_sign = sign(zk);
zk_bits = (zk_sign>0);

%It's finding the wrong one (There are two pilots, but a lock on second)
[corr_frame, corr_tau_frame] = xcorr(pilot, zk_bits); %locking in on the wrong one
[~, offset_frame] = max(abs(corr_frame)); %locking in on the wrong one
tau_frame = abs(corr_tau_frame(offset_frame))+1;

msg_start = tau_frame-(length(pilot)/2); %Pushes to be +length(pilot) beyond msg???

pilot_start = tau_frame-(3*length(pilot)/2); %Aligns to Pilot start

msg_eye = zk(msg_start:tau_frame-1); %Message found and located

pilot_eye = zk(pilot_start:tau_frame-(length(pilot)/2)-1); %Should line-up to pilot end before msg

zk_bits2 = zk_bits(msg_start:end);

[corr_frame_end, corr_tau_frame_end] = xcorr(pilot, zk_bits2);
[~, offset_end] = max(abs(corr_frame_end));
tau_frame_end = abs(corr_tau_frame_end(offset_end)); %Determines end of msg

zk_msg = zk(msg_start:msg_start+tau_frame_end-1); %comes out to zk(41:50)

ho_hat = (pilot * pilot_eye') / (norm(pilot)^2);

%% Detect bits - One Tap Channel

%zk_msg = zk(41:50);

vk = zk_msg / ho_hat;

xk_hat = sign(vk);
xk_hat = (xk_hat>0);
xk_hat = reshape(xk_hat,1,length(xk_hat));

vIk = vk(1:2:end);
vQk = vk(2:2:end);

v_k = vIk + j*vQk;

% Compute Bit error rate (BER)
BER = mean(xk_hat ~= msg);
disp(['BER = ' num2str(BER)])
disp(' ')

%% Define Constellation
qam_range = 1:sqrt(qam);
qam_range = d*qam_range - 0.5*d - sqrt(qam)/2;
constellation = [];
 
for xi = qam_range
    for xq = qam_range
        constellation = [constellation, xi + 1i*xq];
    end
end


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
    plot(y_received_timing)
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
    stem(1:length(vk),vk,'b')
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

% one of these should provide an AWGN channel coming back to us which is a
% good thing.