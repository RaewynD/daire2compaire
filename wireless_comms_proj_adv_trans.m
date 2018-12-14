%% Raewyn Duvall and Emmanuel Aire-Oaihimire 
%  Team: Daire2Compaire
%  18-758 Wireless Communications
%  Fall 2018


%% --Main Transmit Code-- %%

clear
close all
%clc
% Set seed for random number generator (for repeatability of simulation)
rng('default');

global freq_preamble timing_preamble pilot_size spreading_gain

%user defined values
picture = 100;
srrc = 1;
showplot = 1;
freq_preamble = 4;
timing_preamble = 30;
pilot_size = 4;
msg_size = 106;
delay_size = 50;
spreading_gain = 100;

pwr = 0.075;

d = 1;

fs = 200e6; %sampling
Ts = 1/fs;

fc = 11.25e6; %carrier= 11.25MHz therefore Nyquist is 22.5MHz, which is below the hardware limit
Tc = 1/fc;

T_sym = 0.1e-6; %100*Tn; %sec / symbol
F_sym = 1/T_sym;

symLen = T_sym * fs; %samples per symbol

N = 8; %length of filter in symbol periods (0-24)
Ns = T_sym*N*fs; %Number of filter samples

a = 0.2; %sigma

time = 800e-6;
bits_per_sym = 2;
transmit_size = (time / T_sym) * bits_per_sym;


%% Establish Filter
% Use sqrt-raised cosine filter form 
% ww=FIRRCOS(N,Fc,R,Fs,'rolloff',TYPE) as p^b(t)
if srrc == 1
    p = firrcos(Ns, 1/2/T_sym, a, fs, 'rolloff', 'sqrt');
    % '1/fs' simply serves as 'delta' to approximate integral as sum
    p = p/norm(p)/sqrt(1/fs);
% Use rectangular pulse as filter
else
%    p = [zeros(ceil(Ns/2-T_sym*fs/2),1);
%         ones(ceil(T_sym*fs),1);
%         zeros(N-T_sym*fs-ceil(Ns/2-T_sym*fs/2),1)];
%    p = p/norm(p)/sqrt(1/fs);
    
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
[bits,imdim] = get_pic(picture);
imrecon = reshape(bits,imdim);

msg = bits;
len = length(msg);

% make sure msg is factor of message_size
to_add = zeros(1, msg_size - mod(len, msg_size));
msg = [msg,to_add];

num_msg = length(msg)/msg_size;

x1 = [freq, timing, pilot];

timing_neg = 2*timing-1;
pilot_neg = 2*pilot-1;
pilot_plot = [freq, timing_neg, pilot_neg];

for x = [0:num_msg-1]
    mini_msg = msg((msg_size*x)+1:msg_size*(x+1));
    x1 = [x1, mini_msg, pilot];
    pilot_plot = [pilot_plot, zeros(1,msg_size), pilot_neg];
end

len = length(x1);
if (len > transmit_size)
    disp('Your transmit signal is too long.');
end

% make 1 and -1
x2 = 2*x1-1;
x2 = x2';

pilot_plot = pilot_plot';

% make scaled to qam
xI_base = x2(1:2:end);
xQ_base = x2(2:2:end);

x_base = xI_base + j * xQ_base;
[x_spread, spreading_mask] = spread_msg(x_base);

freq_spread_start = 1;
freq_spread_end = (freq_preamble/2)*spreading_gain - 1;
freq_spread = x_spread(freq_spread_start : freq_spread_end);

timing_spread_start = freq_spread_end + 1;
timing_spread_end = timing_spread_start + (timing_preamble/2)*spreading_gain - 1;
timing_spread = x_spread(timing_spread_start : timing_spread_end);

pilot_spread_start = timing_spread_end + 1;
pilot_spread_end = pilot_spread_start + (pilot_size/2)*spreading_gain - 1;
pilot_spread = x_spread(pilot_spread_start : pilot_spread_end);

msg_spread = zeros((msg_size/2)*spreading_gain,1);

pilot_2_start = pilot_spread_end + (msg_size/2)*spreading_gain + 1;
pilot_2_end = pilot_2_start + (pilot_size/2)*spreading_gain - 1;
pilot_2 = x_spread(pilot_2_start : pilot_2_end);

pilot_plot = pilot_plot(1:2:end) + j * pilot_plot(2:2:end);
pilot_plot = upsample(pilot_plot,spreading_gain);

% convolve with pulse
x_up = upsample(x_spread, fs/F_sym);
x = conv(x_up, p);

% TODO - remove
pilot_plot = [freq_spread;timing_spread;pilot_spread;msg_spread;pilot_2];
pilot_plot = upsample(pilot_plot, fs/F_sym);
pilot_plot = conv(pilot_plot,p);

%% transmit complex symbols
transmitsignal = x;
transmitsignal = transmitsignal/max(abs(transmitsignal));

rand2 = ceil(rand([1,1])*delay_size/2)*2 + delay_size;
rand3 = ceil(rand([1,1])*delay_size/2)*2 + delay_size*3;
rand4 = ceil(rand([1,1])*delay_size/2)*2 + (delay_size*5);

transmitsignal1 = [transmitsignal; zeros(rand4,1)]*pwr;

transmitsignal2 = [zeros(rand2,1); transmitsignal; zeros(rand4-rand2,1)]*pwr*0.95;

transmitsignal3 = [zeros(rand3,1); transmitsignal; zeros(rand4-rand3,1)]*pwr*0.9;

transmitsignal4 = [zeros(rand4,1); transmitsignal]*pwr*0.85;

%transmitsignal = (transmitsignal1 + transmitsignal2 + transmitsignal3)/3;% + transmitsignal4)/4;

transmitsignal = reshape(transmitsignal, [], 1);

pilot_plot = pilot_plot/max(abs(pilot_plot));
pilot_plot = pilot_plot*pwr;
pilot_plot = reshape(pilot_plot, [], 1);

freq_spread_len = length(freq_spread);
timing_spread_len = length(timing_spread);
pilot_spread_len = length(pilot_spread);
msg_spread_len = length(msg_spread);

save('transmitsignal.mat','transmitsignal')

% save for analysis in receive
save global_vars.mat d fs Ts fc Tc T_sym F_sym symLen a p timing pilot msg ...
    N Ns num_msg pilot_plot bits imdim msg_size ...
    spreading_gain spreading_mask freq_spread_len timing_spread ...
    timing_spread_len pilot_spread_len msg_spread_len

if srrc == 1
    save('transmitsignal_SRRC.mat','transmitsignal')
elseif srrc == 0
    save('transmitsignal_RECT.mat','transmitsignal')
else
end

%% Plot time and frequency domain signals
%load receivedsignal.mat

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
    h(1) = subplot(3,2,1);
    plot([1:length(p)]/fs,p)
    ylabel('$p^{transmit}(t)$')
    xlabel('Time (s)')
    title('Pulse Signal')
    set(gca,'fontsize', 15)
    ax(1) = subplot(3,2,3);
    h(3) = ax(1);
    plot(real(transmitsignal),'b')
    hold on
    plot(imag(transmitsignal),'r')
    plot(real(pilot_plot),'g')
    plot(imag(pilot_plot),'y--')
    legend('real','imag')
    ylabel('$x^{I}(t)$,  $x^{Q}(t)$')
    xlabel('Time in samples')
    title('Transmitted Signal')
    set(gca,'fontsize', 15)
    %%% please add this as a plot
    
    figure(3)
    plot(real(transmitsignal1),'Color',[0,0,0.7])
    hold on
    plot(imag(transmitsignal1),'Color',[0,0,0.5])
    plot(real(transmitsignal2),'Color',[0,0.7,0])
    plot(imag(transmitsignal2),'Color',[0,0.5,0])
    plot(real(transmitsignal3),'Color',[0.7,0,0])
    plot(imag(transmitsignal3),'Color',[0.5,0,0])
    plot(real(transmitsignal4),'Color',[0,0.7,0.7])
    plot(imag(transmitsignal4),'Color',[0,0.5,0.5])
    
    figure(1)
    hold on
    h(2) = subplot(3,2,2);
    plot([-lenp/2+1:lenp/2]/lenp*fs,20*log10(abs(fftshift(1/sqrt(lenp)*fft(p)))))
    ylabel('$|P^{transmit}(f)|$')
    xlabel('Frequency')
    title('Frequency Response of Pulse')
    axis([-4*fc 4*fc -Inf Inf])
    set(gca,'fontsize', 15)
    ax(2) = subplot(3,2,4);
    h(4) = ax(2);
    plot([0:length(transmitsignal)-1]/length(transmitsignal)-0.5, abs(fftshift(fft(transmitsignal))))
    ylabel('$|X^{base}(f)|$')
    xlabel('Frequency in 1/samples')
    title('Frequency Response of Transmitted Signal')
    set(gca,'fontsize', 15)
    %linkaxes(ax,'x')
    zoom on
    h(5) = subplot(3,2,5);
    imshow(imrecon)
    xlabel('Desired Image')
    pos = get(h,'Position');
    new = mean(cellfun(@(v)v(1),pos(1:2)));
    set(h(5),'Position',[new,pos{end}(2:end)])
    set(gca,'fontsize', 15)
%     bits_str = mat2str(bits);
%     msg_str = mat2str(msg);
%     %str = ['bits: \n', bits_str, '\n', 'msg: \n', msg_str];
%     ax = subplot(3,2,6);
%     txt = 1;
%     text(0,txt,'bits');
%     txt = txt - 0.2;
%     bits_txt = reshape(bits,imdim);
%     for x = 1:imdim(1)
%         str = mat2str(bits_txt(x,:));
%         text(0,txt,str);
%         txt = txt - 0.2;
%     end
%     text(0,txt,'coded msg');
%     txt = txt - 0.2;
%     msg_txt = reshape(test,[imdim(1),imdim(2)*2]);
%     for x = 1:imdim(1)
%         str = mat2str(msg_txt(x,:));
%         text(0,txt,str);
%         txt = txt - 0.2;
%     end
%     set(ax, 'visible', 'off')
    
    
else
    close all
end

%% ---Helper Functions--- %%

% spread code
function [spreaded_sym,spread_mask] = spread_msg(msg)

    global spreading_gain

    len = length(msg);
    
    spreaded_sym = [];
    spread_mask = [];
    
    for x = 1:len
        spread = get_bits(3)';
        spread = spread * 2 - 1;
        
        spread = msg(x)*spread;
        spreaded_sym = [spreaded_sym;spread];
        spread_mask = [spread_mask;spread];
    end

end

% get frequency preamble, timing preamble, and pilot to transmit
function bits = get_bits(pic)

    global freq_preamble timing_preamble pilot_size spreading_gain

    switch pic
        case 0
            bits = ones(1,freq_preamble);
        case 1
            bits = randi([0 1],1,timing_preamble);
        case 2
            bits = randi([0 1],1,pilot_size);
        case 3
            bits = randi([0 1],1,spreading_gain);
        otherwise
            bits = [1 1 0 1 0 0 0 1 0 1];
    end
end

function [bits,imdim] = get_pic(pic)

    switch pic
        case 4
            A = [1,1;0,1];
            imdim = size(A);
            bits = A(:);
            bits = bits';
        case 88
            A = imread('shannon88.bmp');
            imdim = size(A);
            bits = A(:);
            bits = bits';
        case 100
            A = [
            1,1,1,0,1,0,1,0,0,0;
            1,0,1,0,1,1,1,0,0,0;
            1,1,1,0,1,0,1,0,0,0;
            0,0,0,1,0,1,0,1,1,1;
            0,0,0,1,1,1,0,0,1,0;
            0,0,0,1,0,1,0,1,1,1;
            1,0,1,0,0,1,0,0,1,0;
            1,0,1,0,1,0,1,0,1,0;
            0,1,0,0,1,1,1,0,0,0;
            0,1,0,0,1,0,1,0,1,0    
            ];
            imdim = size(A);
            bits = A(:);
            bits = bits';
        case 816
            A = imread('shannon816.bmp');
            imdim = size(A);
            bits = A(:);
            bits = bits';
        case 3036
            A = imread('shannon3036.bmp');
            imdim = size(A);
            bits = A(:);
            bits = bits';
        case 6596
            A = imread('shannon6596.bmp');
            imdim = size(A);
            bits = A(:);
            bits = bits';
        case 13720
            A = imread('shannon13720.bmp');
            imdim = size(A);
            bits = A(:);
            bits = bits';
        case 24180
            A = imread('shannon24180.bmp');
            imdim = size(A);
            bits = A(:);
            bits = bits';
        case 46260
            A = imread('shannon46260.bmp');
            imdim = size(A);
            bits = A(:);
            bits = bits';
        otherwise
            bits = [1 1 0 1 0 0 0 1 0 1];
    end
end

