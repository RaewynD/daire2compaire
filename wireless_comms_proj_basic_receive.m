%% Emmanuel Aire-Oaihimire and Raewyn Duvall
%  Team: Daire2Compaire
%  18-758 Wireless Communications
%  Fall 2018
 
%% --Main Receive Code-- %%

clear
close all
clc

% Define User Values
srrc = 0;
real_time = 1;

if srrc == 1
    load receivedsignal_SRRC
    load transmitsignal_SRRC
else
    load receivedsignal_RECT
    load transmitsignal_RECT
end

if real_time == 1
    load receivedsignal
end

load global_vars
%d fs Ts fn Tn T_sym F_sym symLen a p timing pilot msg

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

[corr, corr_tau] = xcorr(timing_sent, y_received);
[~, offset] = max(abs(corr));
tau = abs(corr_tau(offset)+1);

y_received = y_received(tau:end);

%% Grab and separate into REAL and IMAGINARY

yI = real(y_received);
yQ = imag(y_received);

%% Filter low pass signals with matched filter in each arm

% '1/fs' simply serves as 'delta' to approximate integral as sum
zI = conv(w,yI)*(1/fs);
zQ = conv(w,yQ)*(1/fs);

%% Sample filtered signal

k = 1;
zIk = [];
zQk = [];
for s = length(w)+1:T_sym/Ts:length(zI)
    zIk(k) = zI(s);
    zQk(k) = zQ(s);
    k = k+1;
end

%% Frame Recovery

%{
fIk = frame(1:2:end);
fQk = frame(2:2:end);
len = min([length(fIk) length(fQk)]);

zI_bits = sign(zIk); 
zQ_bits = sign(zQk);
zIbits = (zI_bits>0)';
zQbits = (zQ_bits>0)';

max_f = 1;
max_count = 0;
for f = 1:(length(zIbits) - len)
    zI_check = zIbits(f:f+len-1);
    zQ_check = zQbits(f:f+len-1);
    count = 0;
    for l = 1:len
        if (zI_check(l) == fIk(l))
            count = count + 1;
        end
        if (zQ_check(l) == fQk(l))
            count = count + 1;
        end
    end
    if count > max_count
        max_count = count;
        max_f = f;
    end
end
%}


%% Find Message

known_bits = [1 1 0 1 0 0 0 1 0 1];

xIk = known_bits(1:2:end);
xQk = known_bits(2:2:end);
len = min([length(xIk) length(xQk)]);

%zI_bits = sign(zIk); 
%zQ_bits = sign(zQk);
%zIbits = (zI_bits>0)';
%zQbits = (zQ_bits>0)';

k = strfind(zIbits, xIk);
l = strfind(zQbits, xQk);

overlap = ismember(k,l);

bits = 1:(min([length(zIbits) length(zQbits)])*2);
for x = 1:(min([length(zIbits) length(zQbits)]))
    bits(x*2-1) = zIbits(x);
    bits(x*2) = zQbits(x);
end

%m = strfind(bits, known_bits);

%if (isempty(m))
%    disp('Not in full bits...')
%end

zIk_frame = zIk(k:k+len-1);
zQk_frame = zQk(l:l+len-1);

hoI_hat = (xIk * zIk_frame) / (norm(xIk)^2);
hoQ_hat = (xQk * zQk_frame) / (norm(xQk)^2);

%% Detect bits - One Tap Channel

vIk = zIk_frame / hoI_hat;
vQk = zQk_frame / hoQ_hat;

vk = vIk + j*vQk;

xIk_hat = sign(vIk); 
xQk_hat = sign(vQk);
bitI_hat = (xIk_hat>0);
bitQ_hat = (xQk_hat>0);
bits_hat = reshape([bitI_hat'; bitQ_hat'],1,L)

% Compute Bit error rate (BER)
BER = mean(bits_hat ~= known_bits);
disp(['BER = ' num2str(BER)])
disp(' ')

%% Define Constellation
qam_range = 1:sqrt(qam);
qam_range = D*qam_range - 0.5*D - sqrt(qam)/2;
constellation = [];
 
for xi = qam_range
    for xq = qam_range
        constellation = [constellation, xi + 1i*xq];
    end
end

% Part C - make decision
%x_hat = zeros(L,1);

%P_error = 0;
%for x_count = 1:L
%    [min_dist, index] = min(abs(constellation - x_received(x_count)));
%    x_hat(x_count) = constellation(index);
%    if(x_hat(x_count) ~= x_transmitted(x_count))
%        P_error = P_error + 1;
%    end
%end

% Calculate and print error performance
%Pe = P_error / L;
%fprintf('Symbols wrong: %f Symbol error rate Pe = %f\n',P_error, Pe)

if graph == 1
    constellationmarkersize = 6;
    ax = [];
    cats = 0;
    lono = 1:x;

    close all
    figure(1)
    LargeFigure(gcf, 0.15); % Make figure large
    clf
    zoom off
    plot(constellation,'rs','MarkerSize',constellationmarkersize,'MarkerFaceColor','r')
    set(gca,'DataAspectRatio',[1 1 1])
    grid on
    hold on
    D = max(D, max(abs(vk))+1);
    axis([-D D -D D])
    plot([-D:D/100:D],zeros(size([-D:D/100:D])),'k','LineWidth',2)
    plot(zeros(size([-D:D/100:D])),[-D:D/100:D],'k','LineWidth',2)
    set(gca,'fontsize', 15)
    xlabel('$x^{I}$, $z^{I}$')
    ylabel('$x^{Q}$, $z^{Q}$')

    title('Constellation Plot')
    plot(constellation,'rs','MarkerSize',constellationmarkersize,'MarkerFaceColor','r')
    for ii=1:L/2
        plot(vk(ii),'bx')
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
    plot([1:length(p)]/fs,p)
    ylabel('$p^{transmit}(t)$')
    title('Pulse Signal p(t)')
    set(gca,'fontsize', 15)
    subplot(3,2,3);
    plot(real(transmitsignal),'b')
    hold on
    plot(imag(transmitsignal),'r')
    legend('real','imag')
    ylabel('$x^{I}(t),    x^{Q}(t)$')
    xlabel('Time in samples')
    title('Transmit Signal')
    set(gca,'fontsize', 15)
    subplot(3,2,2);
    plot([-lenp/2+1:lenp/2]/lenp*fs,20*log10(abs(fftshift(1/sqrt(lenp)*fft(p)))))
    ylabel('$|P^{transmit}(f)|$')
    title('Pulse Signal in Frequency Domain')
    axis([-4*fc 4*fc -40 40])
    set(gca,'fontsize', 15)
    subplot(3,2,4);
    plot([0:length(transmitsignal)-1]/length(transmitsignal)-0.5, abs(fftshift(fft(transmitsignal))))
    ylabel('$|X^{base}(f)|$')
    xlabel('Frequency in 1/samples')
    title('Transimt Signal in Frequency Domain')
    set(gca,'fontsize', 15)
    subplot(3,2,5)
    plot(real(receivedsignal),'b')
    hold on
    plot(imag(receivedsignal),'r')
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
    
    figure(3)
    LargeFigure(gcf, 0.15); % Make figure large
    clf
    ax(1) = subplot(2,2,1);
    %stem([1:x],bitI_hat,'b')
    hold on
    stem(lono,zI_bits,'r')
    for zee = 1:length(zIk_frame)
        cat = find(zIk==zIk_frame(zee));
        stem(lono(cat),zI_bits(cat),'b')
    end
    ylabel('$z^I_{k}$') %'$x^I_k,   z^I_{k}$'
    xlabel('discrete time  $k$  (sampled at $t=kT$)')
    title('Sampler Output $z^I_{k}$')
    ylim([-2 2]);
    set(gca,'fontsize', 15)
    ax(2) = subplot(2,2,2);
    %stem([1:x],bitQ_hat,'b')
    hold on
    stem(lono,zQ_bits,'b')
    for zee = 1:length(zQk_frame)
        cat = find(zQk==zQk_frame(zee));
        stem(lono(cat),zQ_bits(cat),'r')
    end
    ylabel('$z^Q_{k}$') %'$x^Q_k,   z^Q_{k}$'
    xlabel('discrete time  $k$  (sampled at $t=kT$)')
    title('Sampler Output $z^Q_{k}$')
    ylim([-2 2]);
    set(gca,'fontsize', 15)
    subplot(2,2,3);
    stem([1:len],xIk_hat','b')
    hold on
    stem([1:len],xQk_hat','r')
    ylabel('$v^I_{k}, v^Q_{k}$') %'$v^I_k,   v^I_{k}$'
    xlabel('discrete time  $k$  (sampled at $t=kT$)')
    title('Equalizer $v_{k}$ output samples')
    set(gca,'fontsize', 15)
    ylim([-2 2]);
    subplot(2,2,4);
    %stem([1:x],bitQ_hat,'b')
    %hold on
    stem([1:L],bits_hat','k')
    ylabel('$z_{k}$') %'$x^Q_k,   z^Q_{k}$'
    xlabel('discrete time  $k$  (sampled at $t=kT$)')
    title('Decoded Bits $z_{k}$')
    ylim([-2 2]);
    set(gca,'fontsize', 15)
    linkaxes(ax,'x')
    zoom xon

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