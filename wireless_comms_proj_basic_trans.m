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

transmitsignal = (xI + j*xQ);
transmitsignal = reshape(transmitsignal, [], 1);

save('transmitsignal.mat','transmitsignal')

if (srrc == 1)
    save('transmitsignal_SRRC.mat','transmitsignal')
else
    save('transmitsignal_RECT.mat','transmitsignal')
end

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
