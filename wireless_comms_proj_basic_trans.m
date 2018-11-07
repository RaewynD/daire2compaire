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
srrc = 0;

d = 1;

fs = 200e6; %sampling
Ts = 1/fs;

fn = 12.5e6; %nyquist
Tn = 1/fn;

T_sym = 100*Tn; %sec / symbol
F_sym = 1/T_sym;

symLen = T_sym * fs; %samples per symbol

a = 0.2;

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

%% Define binary transmission
freq = get_bits(0);
timing = get_bits(1);
pilot = get_bits(1);
msg = get_bits(picture);
% make sure even number of bits so real and img are equal length
if mod(length(msg), 2) ~= 0
    msg = [msg, 0];
end

x1 = [freq, timing, pilot, msg, pilot];

% make 1 and -1
x2 = 2*x1-1;
x2 = x2';

% make scaled to qam
xI_base = x2(1:2:end)*(0.5*d);
xQ_base = x2(2:2:end)*(0.5*d);

% convolve with pulse
xI_up = upsample(xI_base, fs/F_sym);
xQ_up = upsample(xQ_base, fs/F_sym);
xI = conv(xI_up, p);
xQ = conv(xQ_up, p);

%% transmit complex symbols
transmitsignal = (xI + j*xQ);
transmitsignal = reshape(transmitsignal, [], 1);

save('transmitsignal.mat','transmitsignal')

% save for analysis in receive
save transmitpreamble.mat timing pilot msg

if srrc == 1
    save('transmitsignal_SRRC.mat','transmitsignal')
else
    save('transmitsignal_RECT.mat','transmitsignal')
end

%% Plot time and frequency domain signals

if srrc == 1
    load receivedsignal_SRRC
    load transmitsignal_SRRC
else
    load receivedsignal_RECT
    load transmitsignal_RECT
end

figure(1)
LargeFigure(gcf, 0.15); % Make figure large
clf
subplot(3,2,1);
plot([1:length(p)]/fs,p)
ylabel('$p^{transmit}(t)$')
set(gca,'fontsize', 15)
subplot(3,2,3);
plot(real(transmitsignal),'b')
hold on
plot(imag(transmitsignal),'r')
legend('real','imag')
ylabel('$x^{I}(t)$,  $x^{Q}(t)$')
xlabel('Time in samples')
set(gca,'fontsize', 15)
subplot(3,2,2);
plot([-lenp/2+1:lenp/2]/lenp*fs,20*log10(abs(fftshift(1/sqrt(lenp)*fft(p)))))
ylabel('$|P^{transmit}(f)|$')
axis([-4*fc 4*fc -40 40])
set(gca,'fontsize', 15)
subplot(3,2,4);
plot([0:length(transmitsignal)-1]/length(transmitsignal)-0.5, abs(fftshift(fft(transmitsignal))))
ylabel('$|X^{base}(f)|$')
xlabel('Frequency in 1/samples')
set(gca,'fontsize', 15)
subplot(3,2,5)
plot(real(receivedsignal),'b')
hold on
plot(imag(receivedsignal),'r')
zoom xon
legend('real','imag')
ylabel('$y^{I}(t)$,  $y^{Q}(t)$')
xlabel('Time in samples')
set(gca,'fontsize', 15)
subplot(3,2,6)
plot([0:length(receivedsignal)-1]/length(receivedsignal)-0.5, abs(fftshift(fft(receivedsignal))))
ylabel('$|Y^{base}(f)|$')
xlabel('Frequency in 1/samples')
set(gca,'fontsize', 15)
%subplot(4,2,5);
%stem([1:length(xI_up)],xI_up,'b')
%hold on
%stem([1:length(xI_up)],xQ_up,'r')
%ylabel('$x^I_k,   x^Q_{k}$')
%xlabel('DT transmisison $k$  (sampled at $t=kT$) [before convolution with p(t)]')
%set(gca,'fontsize', 15)
%zoom xon


%% ---Helper Functions--- %%

% get bit array from picture to transmit
function bits = get_bits(pic)
    switch pic
        case 0
            bits = ones(1,50);
        case 1
            bits = randi([0 1],1,50);
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


