%% Raewyn Duvall and Emmanuel Aire-Oaihimire 
%  Team: Daire2Compaire
%  18-758 Wireless Communications
%  Fall 2018


%% --Main Transmit Code-- %%

clear
%close all
%clc
% Set seed for random number generator (for repeatability of simulation)
rng('default');

global freq_preamble timing_preamble pilot_size

%user defined values
picture = 13720;
srrc = 1;
showplot = 1;
freq_preamble = 300;
timing_preamble = 50;
pilot_size = 20;
msg_size = 150;

d = 1;

fs = 200e6; %sampling
Ts = 1/fs;

fc = 11.25e6; %carrier= 11.25MHz therefore Nyquist is 22.5MHz, which is below the hardware limit
Tc = 1/fc;

T_sym = 0.1e-6; %100*Tn; %sec / symbol
F_sym = 1/T_sym;

symLen = T_sym * fs; %samples per symbol

N = 10; %length of filter in symbol periods (0-24
Ns = T_sym*N*fs; %Number of filter samples

a = 0.2; %sigma

time = 800e-6;
bits_per_sym = 2;
transmit_size = floor((time / T_sym) * bits_per_sym)-10;


%% Establish Filter
% Use sqrt-raised cosine filter form 
% ww=FIRRCOS(N,Fc,R,Fs,'rolloff',TYPE) as p^b(t)
if srrc == 1
    p = firrcos(Ns, 1/2/T_sym, a, fs, 'rolloff', 'sqrt');
    % '1/fs' simply serves as 'delta' to approximate integral as sum
    p = p/norm(p)/sqrt(1/fs);
% Use rectangular pulse as filter
else
    p = [zeros(ceil(symLen/2),1);
         ones(ceil(symLen),1);
         zeros(ceil(symLen/2),1)];
    p = p/norm(p)/sqrt(1/fs);
end
% As it curently stands, the length of the pulse is at the value of T,
% which should be Nyquist.

% should be same as symLen
lenp = length(p);

%figure(1)
%plot([1:lenp]/fs,p)

%% Define binary transmission
freq = get_bits(0);
timing = get_bits(1);
pilot = get_bits(2);
msg = get_bits(picture);
len = length(msg);

% make sure msg is factor of message_size
to_add = zeros(1,mod(len, msg_size));
msg = [msg,to_add];

num_msg = length(msg)/msg_size;

x1 = [freq, timing, pilot];

timing_neg = 2*timing-1;
pilot_neg = 2*pilot-1;
pilot_plot = [freq, timing_neg, pilot_neg];

for x = [0:num_msg-1]
    x1 = [x1, msg((msg_size*x)+1:msg_size*(x+1)), pilot];
    pilot_plot = [pilot_plot, zeros(1,msg_size), pilot_neg];
end

len = length(x1)
if (len >= transmit_size)
    disp("Your transmit signal is too long.");
else
    x1 = [ones(1,transmit_size-len),x1];
    pilot_plot = [ones(1,transmit_size-len),pilot_plot];
end
len = length(x1)

% make 1 and -1
x2 = 2*x1-1;
x2 = x2';

pilot_plot = pilot_plot';

% make scaled to qam
xI_base = x2(1:2:end);
xQ_base = x2(2:2:end);

pilot_plot_I = pilot_plot(1:2:end);
pilot_plot_Q = pilot_plot(2:2:end);

% convolve with pulse
xI_up = upsample(xI_base, fs/F_sym);
xQ_up = upsample(xQ_base, fs/F_sym);
xI = conv(xI_up, p);
xQ = conv(xQ_up, p);

pilot_plot_I = upsample(pilot_plot_I, fs/F_sym);
pilot_plot_Q = upsample(pilot_plot_Q, fs/F_sym);
pilot_plot_I = conv(pilot_plot_I, p);
pilot_plot_Q = conv(pilot_plot_Q, p);

%% transmit complex symbols
transmitsignal = (xI + j*xQ);
transmitsignal = transmitsignal/max(abs(transmitsignal));
transmitsignal = reshape(transmitsignal, [], 1);

pilot_plot = (pilot_plot_I + j*pilot_plot_Q);
pilot_plot = pilot_plot/max(abs(pilot_plot));
pilot_plot = reshape(pilot_plot, [], 1);

save('transmitsignal.mat','transmitsignal')

% save for analysis in receive
save global_vars.mat d fs Ts fc Tc T_sym F_sym symLen a p timing pilot msg

if srrc == 1
    save('transmitsignal_SRRC.mat','transmitsignal')
elseif srrc == 0
    save('transmitsignal_RECT.mat','transmitsignal')
else
end

%% Plot time and frequency domain signals
load receivedsignal.mat

ax = []; %Axes connections

if srrc == 1
    %load receivedsignal_SRRC
    load transmitsignal_SRRC
elseif srrc == 0
    %load receivedsignal_RECT
    load transmitsignal_RECT
else
end

if showplot == 1
    figure(1)
    LargeFigure(gcf, 0.15); % Make figure large
    clf
    subplot(2,2,1);
    plot([1:length(p)]/fs,p)
    ylabel('$p^{transmit}(t)$')
    xlabel('Time (s)')
    title('Pulse Signal')
    set(gca,'fontsize', 15)
    ax(1) = subplot(2,2,3);
    plot(real(transmitsignal),'b')
    hold on
    plot(imag(transmitsignal),'r')
    plot(imag(pilot_plot),'y')
    plot(real(pilot_plot),'g')
    legend('real','imag')
    ylabel('$x^{I}(t)$,  $x^{Q}(t)$')
    xlabel('Time in samples')
    title('Transmitted Signal')
    set(gca,'fontsize', 15)
    subplot(2,2,2);
    plot([-lenp/2+1:lenp/2]/lenp*fs,20*log10(abs(fftshift(1/sqrt(lenp)*fft(p)))))
    ylabel('$|P^{transmit}(f)|$')
    xlabel('Frequency')
    title('Frequency Response of Pulse')
    axis([-4*fc 4*fc -Inf Inf])
    set(gca,'fontsize', 15)
    ax(2) = subplot(2,2,4);
    plot([0:length(transmitsignal)-1]/length(transmitsignal)-0.5, abs(fftshift(fft(transmitsignal))))
    ylabel('$|X^{base}(f)|$')
    xlabel('Frequency in 1/samples')
    title('Frequency Response of Transmitted Signal')
    set(gca,'fontsize', 15)
    %linkaxes(ax,'x')
    
    zoom on
else
    close all
end

%% ---Helper Functions--- %%

% get bit array from picture to transmit
function bits = get_bits(pic)

    global freq_preamble timing_preamble pilot_size

    switch pic
        case 0
            bits = ones(1,freq_preamble);
        case 1
            bits = randi([0 1],1,timing_preamble);
        case 2
            bits = randi([0 1],1,pilot_size);
        case 10
            bits = [1 1 0 1 0 0 0 1 0 1];
        case 88
            A = imread('shannon88.bmp');
            bits = A(:);
            bits = bits';
        case 816
            A = imread('shannon816.bmp');
            bits = A(:);
            bits = bits';
        case 3036
            A = imread('shannon3036.bmp');
            bits = A(:);
            bits = bits';
        case 6596
            A = imread('shannon6596.bmp');
            bits = A(:);
            bits = bits';
        case 13720
            A = imread('shannon13720.bmp');
            bits = A(:);
            bits = bits';
        case 24180
            A = imread('shannon24180.bmp');
            bits = A(:);
            bits = bits';
        case 46260
            A = imread('shannon46260.bmp');
            bits = A(:);
            bits = bits';
        otherwise
            bits = [1 1 0 1 0 0 0 1 0 1];
    end
end


