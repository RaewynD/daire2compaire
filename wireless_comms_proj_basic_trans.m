%% --Main Transmit Code-- %%
rng('default'); % Set seed for random number generator (for repeatability of simulation)

%user defined values
picture = 0;
d = 1;

% Define binary transmission
x1 = get_bits(picture);
x_size = size(x1);
x_size = x_size(2);
x_range = 1:ceil(x_size/2);

x2 = x1;
x_transmitted = [];
for x = 1:x_size
    if (x2(x) == 0)
        x2(x) = -1;
    end
end

% translate 2 point pair to 4 qam
for x = x_range
    if (x*2 > x_size)
        x2 = [x2,0];
    end
    x_transmitted = [x_transmitted, x2(x*2 - 1) + (j * x2(x*2))];
end

% give 4 qam proper energy
x_transmitted = 0.5 * d * x_transmitted;

transmitsignal = x_transmitted;

save('transmitsignal.mat','transmitsignal')


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
