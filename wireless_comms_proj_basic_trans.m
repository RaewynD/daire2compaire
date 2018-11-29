%% Raewyn Duvall and Emmanuel Aire-Oaihimire 
%  Team: Daire2Compaire
%  18-758 Wireless Communications
%  Fall 2018


%% --Main Transmit Code-- %%

clear
close all
clc
% Set seed for random number generator (for repeatability of simulation)
rng('default');

%user defined values
picture = -1;
srrc = 2;
showplot = 1;

% MUST BE EVEN
global preamble_size;
preamble_size = 25;

d = 1;

fs = 200e6; %sampling
Ts = 1/fs;

fc = 11.25e6; %carrier= 11.25MHz therefore Nyquist is 22.5MHz, which is below the hardware limit
Tc = 1/fc;

T_sym = 1e-6; %100*Tn; %sec / symbol
F_sym = 1/T_sym;

symLen = T_sym * fs; %samples per symbol

a = 0.2; %sigma

% Use sqrt-raised cosine filter form 
% ww=FIRRCOS(N,Fc,R,Fs,'rolloff',TYPE) as p^b(t)
if srrc == 1
    p = firrcos(Ns, 1/2/T_sym, a, fs, 'rolloff', 'sqrt');
    % '1/fs' simply serves as 'delta' to approximate integral as sum
    p = p/norm(p)/sqrt(1/fs);
% Use rectangular pulse as filter
else
    p = [zeros(ceil(symLen/4),1);
         ones(ceil(symLen/2),1);
         zeros(ceil(symLen/4),1)];
    p = p/norm(p)/sqrt(1/fs);
end

% should be same as symLen
lenp = length(p);

%figure(1)
%plot([1:lenp]/fs,p)

%% Define binary transmission
freq = get_bits(0);
timing = get_bits(1);
pilot = get_bits(1);
msg = get_bits(picture);
% make sure even number of bits so real and img are equal length
if mod(length(msg), 2) ~= 0
    msg = [msg, 0];
end
%frequency sync, timing sync, pilot, msg, pilot
x1 = [freq, timing, pilot, msg, pilot];

% make 1 and -1
x2 = 2*x1-1;
x2 = x2';

% make scaled to qam
xI_base = x2(1:2:end);
xQ_base = x2(2:2:end);

% convolve with pulse
xI_up = upsample(xI_base, fs/F_sym);
xQ_up = upsample(xQ_base, fs/F_sym);
xI = conv(xI_up, p);
xQ = conv(xQ_up, p);

%% transmit complex symbols
transmitsignal = (xI + j*xQ);
transmitsignal = transmitsignal/max(abs(transmitsignal));
transmitsignal = reshape(transmitsignal, [], 1);

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
    load receivedsignal_SRRC
    load transmitsignal_SRRC
elseif srrc == 0
    load receivedsignal_RECT
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
    %axis([-4*fc 4*fc -40 40])
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



%subplot(4,2,5);
%stem([1:length(xI_up)],xI_up,'b')
%hold on
%stem([1:length(xI_up)],xQ_up,'r')
%ylabel('$x^I_k,   x^Q_{k}$')
%xlabel('DT transmisison $k$  (sampled at $t=kT$) [before convolution with p(t)]')
%set(gca,'fontsize', 15)
%zoom xon

%close all

%% ---Helper Functions--- %%

% get bit array from picture to transmit
function bits = get_bits(pic)

    global preamble_size;

    switch pic
        case 0
            bits = ones(1,100);
        case 1
            bits = randi([0 1],1,20);
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


