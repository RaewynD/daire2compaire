%% Raewyn Duvall and Emmanuel Aire-Oaihimire 
%  Team: Daire2Compaire
%  18-758 Wireless Communications
%  Fall 2018


%% --Main Transmit Code-- %%

clear
close all
clc
rng('default'); % Set seed for random number generator (for repeatability of simulation)

%user defined values
picture = 0;
srrc = 0;

d = 1;
T = 1; % Symbol period in microsec
N = 21; % length of filter in symbol periods
alpha = 0.2; % alpha of sqrt raised cosine filter
fc = 12.5; % Carrier Frequency in MHz
fs = 200; % Sampling frequency in MHz
Ns = N*T*fs; % Number of filter samples
% Use sqrt-raised cosine filter form  ww=FIRRCOS(N,Fc,R,Fs,'rolloff',TYPE) as p^b(t)
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

lenp = length(p);

% Define binary transmission
x1 = get_bits(picture);
x_size = size(x1);
x_size = x_size(2);

if mod(x_size,2) ~= 0
    x1 = [x1,0];
end

x2 = 2*x1-1;
x2 = x2';

xI_base = x2(1:2:end)*(0.5*d);
xQ_base = x2(2:2:end)*(0.5*d);

xI_up = upsample(xI_base, fs);
xQ_up = upsample(xQ_base, fs);
xI = conv(xI_up, p);
xQ = conv(xQ_up, p);
len = min([length(xI) length(xQ)]);

transmitsignal = (xI + j*xQ);
transmitsignal = reshape(transmitsignal, [], 1);

save('transmitsignal.mat','transmitsignal')

if srrc == 1
    save('transmitsignal_SRRC.mat','transmitsignal')
else
    save('transmitsignal_RECT.mat','transmitsignal')
end

if srrc == 1
    load receivedsignal_SRRC
    load transmitsignal_SRRC
else
    load receivedsignal_RECT
    load transmitsignal_RECT
end

%% Plot time and frequency domain signals

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
            bits = [1,1,0,1,0,0,0,1,0,1];
    end
end


