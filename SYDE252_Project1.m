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

% plot audio waveform
dt = 1/fs; % time interval (period)
t = 0:dt:(length(bird)- 1)*dt; % all time interval

figure("Name","Audio Waveform Plot");
plot(t, bird);
xlabel('Seconds'); ylabel('Amplitude');

% resample
new_fs = 16000;
upsampled_bird = resample(bird, new_fs, fs); %

% sound(upsampled_bird, new_fs);
% pause(10.5);

dt = 1/new_fs; % time interval (period)
t = 0:dt:(length(upsampled_bird)- 1)*dt; % all time interval

figure("Name","Audio Waveform Plot upsampled");
plot(t, upsampled_bird);
xlabel('Seconds'); ylabel('Amplitude');


%%% part 2: filtering
MA_window_size = 21;
moving_average_filter = movmean(upsampled_bird,[MA_window_size-1 0]);
% sound(moving_average_filter, 16000);


gaussian_filter = [];
Gauss_window_size = 14;
gauss_buffer = zeros([Gauss_window_size-1 1]);
copy_sample1 = cat(1, gauss_buffer, upsampled_bird);

Gauss_factors = gausswin(Gauss_window_size);

for i=Gauss_window_size:size(copy_sample1)
    sample_window = copy_sample1(i-Gauss_window_size+1:i);
    Gauss_sum = sum(sample_window .* Gauss_factors);
    gaussian_filter(i) = Gauss_sum / Gauss_window_size;
end
gaussian_filter = gaussian_filter(Gauss_window_size:end);
% sound(gaussian_filter, 16000);


median_filter = [];
median_window_size = 2;
median_buffer = zeros([median_window_size-1 1]);
copy_sample2 = cat(1, median_buffer, upsampled_bird);

for i=median_window_size:size(copy_sample2)
    sample_window = copy_sample2(i-median_window_size+1:i);
    median_filter(i) = median(sample_window);
end
median_filter = median_filter(median_window_size:end);
% sound(median_filter, 16000);

yma = MA_filter(upsampled_bird, 21, 41);

% silent = 0;
% silent_sound = [];
% for i=1:size(upsampled_bird)
%     if abs(upsampled_bird(i)) < 0.02
%         silent = silent + 1;
%         silent_sound(silent) = upsampled_bird(i);
%     end
% end
% sample_time = length(bird)./fs;
% silent_time = sample_time * (silent/length(upsampled_bird));
% disp(silent_time);
% sound(silent_sound, 16000)

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
silent_time = sample_time * (silent/length(upsampled_bird));
disp(silent_time);

function y = MA_filter(x, L, plot_num)
    [sample_num, ~] = size(x);
    b = (1/L)*ones(1,L);
    y = filter(b, 1, x);
    figure(plot_num);
    plot((1:sample_num), y); % visualize MA output
    axis([0 sample_num -1 1]);
    xlabel("Sample number");
    ylabel("Amplitude");
    title("Moving average filter output");
end