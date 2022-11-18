%%% part 1: preparations

[bird, fs] = audioread("Birds.wav");

% find sampling rate
disp("Sampling rate: " + fs);

% determine stereo or mono
[num_rows, num_columns] = size(bird);
input_channel = "";
if (num_columns == 1); input_channel = "mono"; else; input_channel = "stereo"; end
disp("Input sound file channel: " + input_channel);

% play sound
%sound(bird, fs);

% write sound to a new file
audiowrite("newBirds.wav", bird, fs);

plot_num = 1;
% plot audio waveform
dt = 1/fs; % time interval (period)
t = 0:dt:(length(bird)- 1)*dt; % all time interval

figure(plot_num);
plot_num = plot_num + 1;
plot(t, bird);
xlabel('Seconds'); ylabel('Amplitude');
title("Audio Waveform Plot");

% resample
new_fs = 16000;
x_mono_16k = resample(bird, new_fs, fs); 

% sound(x_mono_16k, new_fs);
% pause(10.5);

dt = 1/new_fs; % time interval (period)
t = 0:dt:(length(x_mono_16k)- 1)*dt; % all time interval

figure(plot_num);
plot_num = plot_num + 1;
plot(t, x_mono_16k);
xlabel('Seconds'); ylabel('Amplitude');
title("Audio Waveform Plot upsampled");





%%% part 2: filtering

% moving average filter
MA_window_size = 21;
y_ma = MA_filter(x_mono_16k, MA_window_size);
% sound(y_ma, 16000);



% weighted moving average filter
gaussian_filter = [];
Gauss_window_size = 14;
gauss_buffer = zeros([Gauss_window_size-1 1]);
copy_sample1 = cat(1, gauss_buffer, x_mono_16k);

Gauss_factors = gausswin(Gauss_window_size);

for i=Gauss_window_size:size(copy_sample1)
    sample_window = copy_sample1(i-Gauss_window_size+1:i);
    Gauss_sum = sum(sample_window .* Gauss_factors);
    gaussian_filter(i) = Gauss_sum / Gauss_window_size;
end
gaussian_filter = gaussian_filter(Gauss_window_size:end);
y_wa = gaussian_filter;
% sound(gaussian_filter, 16000);



% median filter
median_filter = [];
median_window_size = 2;
median_buffer = zeros([median_window_size-1 1]);
copy_sample2 = cat(1, median_buffer, x_mono_16k);

for i=median_window_size:size(copy_sample2)
    sample_window = copy_sample2(i-median_window_size+1:i);
    median_filter(i) = median(sample_window);
end
median_filter = median_filter(median_window_size:end);
y_med = median_filter;
% sound(median_filter, 16000);



%%% part 2.5: find best window size (tiffany's code)

L_MAX = 150;
ma_noise_input = zeros([L_MAX, 1]);
wa_noise_input = zeros([L_MAX, 1]);
med_noise_input = zeros([L_MAX, 1]);

ma_SNR = zeros([L_MAX, 1]);
wa_SNR = zeros([L_MAX, 1]);
med_SNR = zeros([L_MAX, 1]);
for L = 1:L_MAX
    y_ma = MA_filter(x_mono_16k, L);
    y_wa = WA_filter(x_mono_16k, L);
    y_med = MED_filter(x_mono_16k, L);

    ma_noise_input(L) = sum(abs(x_mono_16k-y_ma))/sum(abs(x_mono_16k));
    wa_noise_input(L) = sum(abs(x_mono_16k-y_wa))/sum(abs(x_mono_16k));
    med_noise_input(L) = sum(abs(x_mono_16k-y_med))/sum(abs(x_mono_16k));

    ma_SNR(L) = 20*log10(sum(y_ma.*y_ma)/sum((x_mono_16k-y_ma).*(x_mono_16k-y_ma)));
    wa_SNR(L) = 20*log10(sum(y_wa.*y_wa)/sum((x_mono_16k-y_wa).*(x_mono_16k-y_wa)));
    med_SNR(L) = 20*log10(sum(y_med.*y_med)/sum((x_mono_16k-y_med).*(x_mono_16k-y_med)));
end

figure(plot_num);
plot_num = plot_num + 1;
plot(1:L_MAX, ma_noise_input);
hold on;
plot(1:L_MAX, wa_noise_input);
plot(1:L_MAX, med_noise_input);
legend('Moving average filter', 'Weighted average filter', 'Median filter','Location','southeast');
xlabel("Window size");
ylabel("Noise to input ratio");
title("Noise to input ratio");
hold off;

figure(plot_num);
plot_num = plot_num + 1;
plot(1:L_MAX, ma_SNR);
hold on;
plot(1:L_MAX, wa_SNR);
plot(1:L_MAX, med_SNR);
legend('Moving average filter', 'Weighted average filter', 'Median filter');
xlabel("Window size");
ylabel("SNR (dB)");
title("SNR ratio");
hold off;


%%% part 3: signal analysis
% detect silent region
silent = 0;
silent_counter = 0;
% min silent length = 100 samples, silent threshold: 0 ~ -0.010
for i=1:size(moving_average_filter)
    % start counting silent length
    if moving_average_filter(i) > -0.010 && moving_average_filter(i) < -0.002
        silent_counter = silent_counter + 1;
    else
        if silent_counter > 100
            silent = silent + silent_counter;
        end
        silent_counter = 0;
    end
end
sample_time = length(bird)./fs;
silent_time = sample_time * (silent/length(x_mono_16k));
disp(silent_time);



% tiffany's functions

function y = MA_filter(x, L)
    [sample_num, ~] = size(x);
    b = (1/L)*ones(1,L);
    y = filter(b, 1, x);
    figure(41);
    plot((1:sample_num), y); % visualize MA output
    axis([0 sample_num -1 1]);
    xlabel("Sample number");
    ylabel("Amplitude");
    title("Moving average filter output");
end

% Weighted average filter
function y = WA_filter(x, L)
    b = gausswin(L);
    b_norm = b / sum(b);
    y = filter(b_norm, 1, x);
    % wvtool(gausswin(L));
end

% Median filter
function y = MED_filter(x, L)
    [sample_num, ~] = size(x); 
    y = zeros(size(x));
    x = cat(1,zeros([L-1, 1]),x);
    for n = 1:sample_num
        y(n) = median(x(n:n+L-1));
    end
end