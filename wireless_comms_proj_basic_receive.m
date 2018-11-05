clear
close all
clc

srrc = 1;

if srrc == 1
    load receivedsignal_SRRC
    load transmitsignal_SRRC
else
    load receivedsignal_RECT
    load transmitsignal_RECT
end

x_received = receivedsignal;
x_transmitted = transmitsignal;

%% Grab and separate into REAL and IMAGINARY

yI = real(x_received);
yQ = imag(x_received);

%% Define User Values
qam = 4;
L = 10; % Number of bits in message
graph = 1;
D = 1;
LL = 200; % Total number of bits
T = 1; % Symbol period in microsec
N = 21; % length of filter in symbol periods
alpha = 0.2; % alpha of sqrt raised cosine filter
fc = 12.5; % Carrier Frequency in MHz
fs = 200; % Sampling frequency in MHz
Ns = N*T*fs; % Number of filter samples

%% Define Pulse

if srrc == 1
    p = firrcos(Ns,1/2/T,alpha,fs,'rolloff','sqrt');
    p = p/norm(p)/sqrt(1/fs); % '1/fs' simply serves as 'delta' to approximate integral as sum
% Use rectangular pulse as one possible filter ???
else
    p = [zeros(ceil(Ns/2-T*fs/2),1);
    ones(T*fs,1);
    zeros(N-T*fs-ceil(Ns/2-T*fs/2),1) ];
    p = p/norm(p)/sqrt(1/fs);
end

%% Define Matched Filter

% Matched filter
w = flipud(p);

% A rectangular (ideal) RF filter of bandwidth 3/T (typically RF filter is quite broadband)
whalflen = 20*fs*T;
wsmoothhalflen = ceil(whalflen/100);


%% Apply Timing Recovery



%% Filter low pass signals with matched filter in each arm
zI = conv(w,yI)*(1/fs); % '1/fs' simply serves as 'delta' to approximate integral as sum
zQ = conv(w,yQ)*(1/fs); % '1/fs' simply serves as 'delta' to approximate integral as sum 


%% Sample filtered signal
zIk = zI(Ns+(2*whalflen)+1:fs*T:end); 
zQk = zQ(Ns+(2*whalflen)+1:fs*T:end); 

%% Frame Recovery

known_bits = [1 1 0 1 0 0 0 1 0 1];

xIk = known_bits(1:2:end);
xQk = known_bits(2:2:end);
len = min([length(xIk) length(xQk)]);

zIbits = sign(zIk); 
zQbits = sign(zQk);
zIbits = (zIbits>0);
zQbits = (zQbits>0);

k = strfind(zIbits', xIk);
l = strfind(zQbits', xQk);

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
bits_hat = reshape([bitI_hat'; bitQ_hat'],L,1)

%% Part A - define constellation
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
    xlabel('x^I, z^I')
    ylabel('x^Q, z^Q')

    title('Constellation Plot')
    plot(constellation,'rs','MarkerSize',constellationmarkersize,'MarkerFaceColor','r')
    for ii=1:L/2
        plot(vk(ii),'bx')
        plot(constellation,'rs','MarkerSize',constellationmarkersize,'MarkerFaceColor','r')
        if (rem(ii,100)==0)
            pause(.00002)
        end
    end


    % Display signals
    figure(2)
    clf
    subplot(2,1,1)
    plot(real(transmitsignal),'b')
    hold on
    plot(imag(transmitsignal),'r')
    set(gca,'fontsize', 15)
    legend('real','imag')
    ylabel('xI(t)  and  xQ(t)')
    xlabel('Time in samples')
    subplot(2,1,2)
    plot(zI,'b')
    hold on
    plot(zQ,'r')
    set(gca,'fontsize', 15)
    zoom xon
    legend('real','imag')
    ylabel('yI(t)  and  yQ(t)')
    xlabel('Time in samples')

    figure(3)
    clf
    subplot(2,1,1)
    plot([0:length(transmitsignal)-1]/length(transmitsignal)-0.5, abs(fftshift(fft(transmitsignal))))
    set(gca,'fontsize', 15)
    ylabel('abs(X(f))')
    xlabel('Frequency in 1/samples')
    subplot(2,1,2)
    plot([0:length(receivedsignal)-1]/length(receivedsignal)-0.5, abs(fftshift(fft(receivedsignal))))
    set(gca,'fontsize', 15)
    ylabel('abs(Y(f))')
    xlabel('Frequency in 1/samples')
    
    figure(4)
    subplot(2,1,1);
    stem([1:L/2],bitI_hat,'b')
    hold on
    stem([1:L/2],zIbits,'r')
    ylabel('$x^I_k,   z^I_{k}$')
    xlabel('discrete time  $k$  (sampled at $t=kT$)')
    subplot(2,1,2);
    stem([1:L/2],bitQ_hat,'b')
    hold on
    stem([1:L/2],zQbits,'r')
    ylabel('$x^Q_k,   z^Q_{k}$')
    xlabel('discrete time  $k$  (sampled at $t=kT$)')
    linkaxes(ax,'x')
    zoom xon

end

%% Things needed on receiving end
% Split into IQ
% Apply Timing Recovery
% Apply Matched Filter
% Apply Sampling for the Zk bits
% Apply some quantizing - Rohit recommends a One-Tap Channel
% Run process to fully lay out the bits

% one of these should provide an AWGN channel coming back to us which is a
% good thing.