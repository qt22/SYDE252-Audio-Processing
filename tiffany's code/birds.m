%%% part 1: preparations

[bird, fs] = audioread("Birds.wav");

% find sampling rate
disp("Sampling rate: " + fs);

% % determine stereo or mono
% [num_rows, num_columns] = size(bird);
% input_channel = "";
% if (num_columns == 1); input_channel = "mono"; else; input_channel = "stereo"; end
% disp("Input sound file channel: " + input_channel);
% 
% % play sound
% %sound(bird, fs);
% 
% % write sound to a new file
% audiowrite("newBirds.wav", bird, fs);
% 
% plot_num = 1;
% % plot audio waveform
% dt = 1/fs; % time interval (period)
% t = 0:dt:(length(bird)- 1)*dt; % all time interval
% 
% figure(plot_num);
% plot_num = plot_num + 1;
% plot(t, bird);
% xlabel('Seconds'); ylabel('Amplitude');
% title("Audio Waveform Plot");
% 
% % resample
% new_fs = 16000;
% x_mono_16k = resample(bird, new_fs, fs); 
% [num_x_16k,~] = size(x_mono_16k);
% 
% % sound(x_mono_16k, new_fs);
% % pause(10.5);
% 
% dt = 1/new_fs; % time interval (period)
% t = 0:dt:(length(x_mono_16k)- 1)*dt; % all time interval
% 
% figure(plot_num);
% plot_num = plot_num + 1;
% plot(t, x_mono_16k);
% xlabel('Seconds'); ylabel('Amplitude');
% title("Audio Waveform Plot upsampled");
% 
% figure(plot_num);
% plot_num = plot_num + 1;
% plot((1:num_rows), bird);
% xlabel('sampling number'); ylabel('Amplitude');
% title("Audio Waveform Plot (sample size")



%%% part 2: filtering


% % Quantitative measures to determine window size
% % Signal - sample number: 6000 - 8000
% % Noise - sample number: 16400 - 17400
% signal_start = 1400;
% signal_end = 5400;
% noise_start = 8000;
% noise_end = 12000;
% x_SNR = 20*log10(sum(bird(signal_start:signal_end).*bird(signal_start:signal_end))/sum(bird(noise_start:noise_end).*bird(noise_start:noise_end)));
% fprintf("Input SNR: %12.12f\n", x_SNR);
% L_MAX = 50; % Ultimate value is 500
% ma_SNR = zeros([L_MAX, 1]);
% wa_SNR = zeros([L_MAX, 1]);
% med_SNR = zeros([L_MAX, 1]);
% for L = 1:L_MAX
%     y_ma = MA_filter(bird, L)*L;
%     y_wa = WA_filter(bird, L)*L;
%     y_med = MED_filter(bird, L)*L;
% 
%     ma_SNR(L) = 20*log10(sum(y_ma(signal_start:signal_end).*y_ma(signal_start:signal_end))/sum(y_ma(noise_start:noise_end).*y_ma(noise_start:noise_end)));
%     wa_SNR(L) = 20*log10(sum(y_wa(signal_start:signal_end).*y_wa(signal_start:signal_end))/sum(y_wa(noise_start:noise_end).*y_wa(noise_start:noise_end)));
%     med_SNR(L) = 20*log10(sum(y_med(signal_start:signal_end).*y_med(signal_start:signal_end))/sum(y_med(noise_start:noise_end).*y_med(noise_start:noise_end)));
% end
% 
% figure(plot_num);
% plot_num = plot_num + 1;
% plot(1:L_MAX, ma_SNR);
% hold on;
% plot(1:L_MAX, wa_SNR);
% plot(1:L_MAX, med_SNR);
% legend('Moving average filter', 'Weighted average filter', 'Median filter');
% xlabel("Window size");
% ylabel("SNR (dB)");
% title("SNR ratio");
% hold off;



% moving average filter
MA_window_size = 21;
y_ma = MA_filter(bird, MA_window_size);
y_ma_best = MA_filter(bird, 47);





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
y_wa = WA_filter(bird, Gauss_window_size);
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





%%% part 3: signal analysis
% detect silent region
silent = 0;
silent_counter = 0;
silent_region = [];
region_counter = 1;
lower_limit = -5.85e-3;
upper_limit = -3.8e-3;

for i=1:size(y_ma)
    % start counting silent length
    if y_ma(i) > lower_limit && y_ma(i) < upper_limit
        silent_counter = silent_counter + 1;
    else
        if silent_counter > 1000
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

% % slient time
% function time = detect_silence(sample, sample_16k, upper_limit, lower_limit, min_size)
%     % min silent length = [min_size] samples, 
%     % silent threshold: lower_limit ~ upper_limit
%     row = 1;
%     record_start = false;
%     start_index = 0;
%     end_index = 0;
% 
%     [sample_size, ~] = size(sample);
% 
%     for i=1:sample_size
%         % start counting silent length
%         if sample(i) > lower_limit && sample(i) < upper_limit
%             if ~record_start
%                 start_index = i;
%                 record_start = true;
%             end
%             end_index = end_index + 1;
%         else
%             if end_index-start_index > min_size
%                 time(row, :) = [start_index end_index];
%                 % disp("silent time intervals (function): " + time(row, :));
%                 row = row + 1;
%                 record_start = false;
%                 start_index = end_index;
%             end
%         end
%     end
%     
%     sample_time = sample_size / 11025;
%     [upsampled_size, ~] = size(sample_16k);
%     silent = 0;
%     
%     silent_time = sample_time * (silent/upsampled_size);
%     disp("silent time (function): " + silent_time);
% 
% 
% end