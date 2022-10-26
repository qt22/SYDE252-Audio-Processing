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
upsampled_bird = resample(bird, 16000, fs);
fs=16000;
sound(upsampled_bird, fs);

dt = 1/fs; % time interval (period)
t = 0:dt:(length(upsampled_bird)- 1)*dt; % all time interval

figure("Name","Audio Waveform Plot upsampled");
plot(t, upsampled_bird);
xlabel('Seconds'); ylabel('Amplitude');