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
[num_x_16k,~] = size(x_mono_16k);

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
y_ma_best = MA_filter(x_mono_16k, 47);
% sound(y_ma, 16000);



% weighted moving average filter
% gaussian_filter = [];
% 
% gauss_buffer = zeros([Gauss_window_size-1 1]);
% copy_sample1 = cat(1, gauss_buffer, x_mono_16k);
% 
% Gauss_factors = gausswin(Gauss_window_size);
% 
% for i=Gauss_window_size:size(copy_sample1)
%     sample_window = copy_sample1(i-Gauss_window_size+1:i);
%     Gauss_sum = sum(sample_window .* Gauss_factors);
%     gaussian_filter(i) = Gauss_sum / Gauss_window_size;
% end
% gaussian_filter = gaussian_filter(Gauss_window_size:end);
Gauss_window_size = 14;
y_wa = WA_filter(x_mono_16k, Gauss_window_size);
y_wa_best = y_wa;
% sound(y_wa, 16000);



% median filter
% median_filter = [];

% median_buffer = zeros([median_window_size-1 1]);
% copy_sample2 = cat(1, median_buffer, x_mono_16k);
% 
% for i=median_window_size:size(copy_sample2)
%     sample_window = copy_sample2(i-median_window_size+1:i);
%     median_filter(i) = median(sample_window);
% end
% median_filter = median_filter(median_window_size:end);
median_window_size = 48;
y_med = MED_filter(x_mono_16k, median_window_size);
y_med_best = y_med;
% sound(y_med, 16000);



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

% Evaluate noise after filtering
% FFT of input - selected region
fft_fs = 8000;
start_samp = 28500;
end_samp = 29500; % num_x_16k;
length = end_samp - start_samp;
freq = fft_fs*(0:floor(length/2))/length; % Get rid of mirrored output
x_fft = fft(x_mono_16k(start_samp:end_samp));
x_fft = mag2db(abs(x_fft));
x_fft = x_fft((1:(floor(length/2))+1));
figure(plot_num);
plot_num = plot_num + 1;
plot(freq, x_fft);
xlabel("Frequency (Hz)");
ylabel("Magnitude (dB)");
title("FFT of audio input");

% FFT of input - entire clip
length = num_x_16k;
x_fft = fft(x_mono_16k);
x_fft = mag2db(abs(x_fft));
x_fft = x_fft((1:(floor(length/2))+1));

% FFT of the output of moving average filter
freq = fft_fs*(0:floor(length/2))/length; % Get rid of mirrored output
ma_fft = fft(y_ma_best);
ma_fft = mag2db(abs(ma_fft));
ma_fft = ma_fft((1:(floor(length/2))+1));
figure(plot_num);
plot_num = plot_num + 1;
plot(freq, ma_fft);
xlabel("Frequency (Hz)");
ylabel("Magnitude (dB)");
title("FFT of moving average output");
figure(plot_num);
plot_num = plot_num + 1;
plot(freq, (ma_fft - x_fft));
xlabel("Frequency (Hz)");
ylabel("Magnitude (dB)");
title("Frequency response of moving average filter");

% FFT of the output of weighted average filter
freq = fft_fs*(0:floor(length/2))/length; % Get rid of mirrored output
wa_fft = fft(y_wa_best);
wa_fft = mag2db(abs(wa_fft));
wa_fft = wa_fft((1:(floor(length/2))+1));
figure(plot_num);
plot_num = plot_num + 1;
plot(freq, wa_fft);
xlabel("Frequency (Hz)");
ylabel("Magnitude (dB)");
title("FFT of weighted average output");
figure(plot_num);
plot_num = plot_num + 1;
plot(freq, (wa_fft - x_fft)); % ERROR: array size too big?!
xlabel("Frequency (Hz)");
ylabel("Magnitude (dB)");
title("Frequency response of weighted average filter");

% FFT of the output of median filter
freq = fft_fs*(0:floor(length/2))/length; % Get rid of mirrored output
med_fft = fft(y_med_best);
med_fft = mag2db(abs(med_fft));
med_fft = med_fft((1:(floor(length/2))+1));
figure(plot_num);
plot_num = plot_num + 1;
plot(freq, med_fft);
xlabel("Frequency (Hz)");
ylabel("Magnitude (dB)");
title("FFT of median output");
figure(plot_num);
plot_num = plot_num + 1;
plot(freq, (med_fft - x_fft)); % ERROR: array size too big?!
xlabel("Frequency (Hz)");
ylabel("Magnitude (dB)");
title("Frequency response of median filter");



%%% part 3: signal analysis
% detect silent region
silent = 0;
silent_counter = 0;
silent_region = [];
region_counter = 1;
lower_limit = -5.8e-3;
upper_limit = -4.1e-3;

for i=1:size(y_ma)
    % start counting silent length
    if y_ma(i) > -5.85e-3 && y_ma(i) < -3.8e-3
        silent_counter = silent_counter + 1;
    else
        if silent_counter > 3000
            silent = silent + silent_counter;
            silent_region(region_counter, :) = [i-silent_counter, i];
            region_counter = region_counter + 1;
        end
        silent_counter = 0;
    end
end
[sample_size, ~] = size(bird);
sample_time = sample_size / fs;
[upsampled_size, ~] = size(x_mono_16k);
silent_time = sample_time * (silent/upsampled_size);
disp("silent time: " + silent_time);

figure(41);
for i=1:size(silent_region)
    silent_interval = silent_region(i, 2) - silent_region(i, 1);
    rectangle('Position', [silent_region(i,1),lower_limit,silent_interval,upper_limit-lower_limit]);
end

% func_silent = detect_silence(y_ma, x_mono_16k, -4.1e-3, -5.8e-3, 1000);



% filter functions

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

% slient time
function time = detect_silence(sample, sample_16k, upper_limit, lower_limit, min_size)
    % min silent length = [min_size] samples, 
    % silent threshold: lower_limit ~ upper_limit
    row = 1;
    record_start = false;
    start_index = 0;
    end_index = 0;

    [sample_size, ~] = size(sample);

    for i=1:sample_size
        % start counting silent length
        if sample(i) > lower_limit && sample(i) < upper_limit
            if ~record_start
                start_index = i;
                record_start = true;
            end
            end_index = end_index + 1;
        else
            if end_index-start_index > min_size
                time(row, :) = [start_index end_index];
                % disp("silent time intervals (function): " + time(row, :));
                row = row + 1;
                record_start = false;
                start_index = end_index;
            end
        end
    end
    
    sample_time = sample_size / 11025;
    [upsampled_size, ~] = size(sample_16k);
    silent = 0;
    
    silent_time = sample_time * (silent/upsampled_size);
    disp("silent time (function): " + silent_time);


end