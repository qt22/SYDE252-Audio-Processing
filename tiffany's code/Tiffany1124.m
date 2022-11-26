% Preparations
plot_num = 1;
[drum, x_drum, x_drum_16k, drum_fs, drum_num_x, ~, ~, plot_num] = preprocessing("Drum.wav", "Drum_out_p1.wav", plot_num);
[bird, x_bird, x_bird_16k, bird_fs, bird_num_x, ~, ~, plot_num] = preprocessing("Birds.wav", "Birds_out_p1.wav", plot_num);
[speech, x_speech, x_speech_16k, speech_fs, speech_num_x, ~, ~, plot_num] = preprocessing("Speech.wav", "Speech_out_p1.wav", plot_num);

% % Quantitative measures to determine window size
% [drum_ma_SNR, drum_wa_SNR, drum_med_SNR, plot_num] = SNR(x_drum, "Drum", 500, 16000, 20000, 43000, 47000, plot_num);
% [bird_ma_SNR, bird_wa_SNR, bird_med_SNR, plot_num] = SNR(x_bird, "Bird", 500, 85000, 87000, 89000, 91000, plot_num);
% [speech_ma_SNR, speech_wa_SNR, speech_med_SNR, plot_num] = SNR(x_speech, "Speech", 500, 15000, 17000, 1, 2001, plot_num);

% Best window size
fprintf("Drum\n");
[~, drum_L_ma_best] = max(drum_ma_SNR);
fprintf("  Moving average filter, best window size: %d\n", drum_L_ma_best);
[~, drum_L_wa_best] = max(drum_wa_SNR);
fprintf("  Weighted average filter, best window size: %d\n", drum_L_wa_best);
[~, drum_L_med_best] = max(drum_med_SNR);
fprintf("Median filter, best window size: %d\n", drum_L_med_best);
fprintf("Bird\n");
[~, bird_L_ma_best] = max(bird_ma_SNR);
fprintf("  Moving average filter, best window size: %d\n", bird_L_ma_best);
[~, bird_L_wa_best] = max(bird_wa_SNR);
fprintf("  Weighted average filter, best window size: %d\n", bird_L_wa_best);
[~, bird_L_med_best] = max(bird_med_SNR);
fprintf("Median filter, best window size: %d\n", bird_L_med_best);
fprintf("Speech\n");
[~, speech_L_ma_best] = max(speech_ma_SNR);
fprintf("  Moving average filter, best window size: %d\n", speech_L_ma_best);
[~, speech_L_wa_best] = max(speech_wa_SNR);
fprintf("  Weighted average filter, best window size: %d\n", speech_L_wa_best);
[~, speech_L_med_best] = max(speech_med_SNR);
fprintf("Median filter, best window size: %d\n", speech_L_med_best);

% Filtering
drum_ma_slight_distort = MA_filter(x_drum, 13);
drum_ma_very_distort = MA_filter(x_drum, 890);
figure(plot_num);
plot_num = plot_num + 1;
plot((1:drum_num_x), x_drum); % visualize input and outputs
axis([0 drum_num_x -1 1]);
xlabel("Sample number");
ylabel("Amplitude");
title("Drum - moving average filter output");
hold on;
plot((1:drum_num_x), drum_ma_slight_distort);
plot((1:drum_num_x), drum_ma_very_distort);
legend('Input', 'Window Size: 13 (Slightly distorted)', 'Window Size: 890 (Very distorted)');
hold off;

bird_ma_slight_distort = MA_filter(x_bird, 7);
bird_ma_very_distort = MA_filter(x_bird, 214);
figure(plot_num);
plot_num = plot_num + 1;
plot((1:bird_num_x), x_bird); % visualize input and outputs
axis([0 bird_num_x -1 1]);
xlabel("Sample number");
ylabel("Amplitude");
title("Bird - moving average filter output");
hold on;
plot((1:bird_num_x), bird_ma_slight_distort);
plot((1:bird_num_x), bird_ma_very_distort);
legend('Input', 'Window Size: 7 (Slightly distorted)', 'Window Size: 214 (Very distorted)');
hold off;

speech_ma_slight_distort = MA_filter(x_speech, 5);
speech_ma_very_distort = MA_filter(x_speech, 1500);
figure(plot_num);
plot_num = plot_num + 1;
plot((1:speech_num_x), x_speech); % visualize input and outputs
axis([0 speech_num_x -1 1]);
xlabel("Sample number");
ylabel("Amplitude");
title("Speech - moving average filter output");
hold on;
plot((1:speech_num_x), speech_ma_slight_distort);
plot((1:speech_num_x), speech_ma_very_distort);
legend('Input', 'Window Size: 5 (Slightly distorted)', 'Window Size: 1500 (Very distorted)');
hold off;

% Weighted average filter
drum_wa_10 = WA_filter(x_drum, 5);
drum_wa_500 = WA_filter(x_drum, 1477);
figure(plot_num);
plot_num = plot_num + 1;
plot((1:drum_num_x), x_drum); % visualize input and outputs
axis([0 drum_num_x -1 1]);
xlabel("Sample number");
ylabel("Amplitude");
title("Drum - weighted average filter output");
hold on;
plot((1:drum_num_x), drum_wa_10);
plot((1:drum_num_x), drum_wa_500);
legend('Input', 'Window Size: 5 (Slightly distorted)', 'Window Size: 1477 (Very distorted)');
hold off;

bird_wa_10 = WA_filter(x_bird, 5);
bird_wa_500 = WA_filter(x_bird, 217);
figure(plot_num);
plot_num = plot_num + 1;
plot((1:bird_num_x), x_bird); % visualize input and outputs
axis([0 bird_num_x -1 1]);
xlabel("Sample number");
ylabel("Amplitude");
title("Bird - weighted average filter output");
hold on;
plot((1:bird_num_x), bird_wa_10);
plot((1:bird_num_x), bird_wa_500);
legend('Input', 'Window Size: 5 (Slightly distorted)', 'Window Size: 217 (Very distorted)');
hold off;

speech_wa_10 = WA_filter(x_speech, 7);
speech_wa_500 = WA_filter(x_speech, 1310);
figure(plot_num);
plot_num = plot_num + 1;
plot((1:speech_num_x), x_speech); % visualize input and outputs
axis([0 speech_num_x -1 1]);
xlabel("Sample number");
ylabel("Amplitude");
title("Speech - weighted average filter output");
hold on;
plot((1:speech_num_x), speech_wa_10);
plot((1:speech_num_x), speech_wa_500);
legend('Input', 'Window Size: 7 (Slightly distorted)', 'Window Size: 1310 (Very distorted)');
hold off;


%TODO: slightly distorted and very distorted window size + plot
%TODO: overlay output to determine the best filter (listen to it)

% Median filter
drum_med_10 = MED_filter(x_drum, 4);
drum_med_500 = MED_filter(x_drum, 43);
figure(plot_num);
plot_num = plot_num + 1;
plot((1:drum_num_x), x_drum); % visualize input and outputs
axis([0 drum_num_x -1 1]);
xlabel("Sample number");
ylabel("Amplitude");
title("Drum - median filter output");
hold on;
plot((1:drum_num_x), drum_med_10);
plot((1:drum_num_x), drum_med_500);
legend('Input', 'Window Size: 4 (Slightly distorted)', 'Window Size: 43 (Very distorted)');
hold off;

bird_med_10 = MA_filter(x_bird, 3);
bird_med_500 = MA_filter(x_bird, 37);
figure(plot_num);
plot_num = plot_num + 1;
plot((1:bird_num_x), x_bird); % visualize input and outputs
axis([0 bird_num_x -1 1]);
xlabel("Sample number");
ylabel("Amplitude");
title("Bird - median filter output");
hold on;
plot((1:bird_num_x), bird_med_10);
plot((1:bird_num_x), bird_med_500);
legend('Input', 'Window Size: 3 (Slightly distorted)', 'Window Size: 37 (Very distorted)');
hold off;

speech_med_10 = MA_filter(x_speech, 3);
speech_med_500 = MA_filter(x_speech, 37);
figure(plot_num);
plot_num = plot_num + 1;
plot((1:speech_num_x), x_speech); % visualize input and outputs
axis([0 speech_num_x -1 1]);
xlabel("Sample number");
ylabel("Amplitude");
title("Speech - median filter output");
hold on;
plot((1:speech_num_x), speech_med_10);
plot((1:speech_num_x), speech_med_500);
legend('Input', 'Window Size: 3 (Slightly distorted)', 'Window Size: 37 (Very distorted)');
hold off;

% % Analyzing signals
% % Drum - detect bpm
% acc1_start = 16000;
% acc1_end = 20000;
% acc2_start = 46000;
% acc2_end = 50000;
% acc3_start = 78000;
% acc3_end = 82000;
% drum_ma_best = MA_filter(x_ma, drum_L_ma_best);
% [~, accent1] = max(drum_ma_best(acc1_start:acc1_end)); % find the peak within the time of interest
% accent1 = accent1 + acc1_start - 1;
% [~, accent2] = max(drum_ma_best(acc2_start:acc2_end));
% accent2 = accent2 + acc2_start - 1;
% [~, accent3] = max(drum_ma_best(acc3_start:acc3_end));
% accent3 = accent3 + acc3_start - 1;
% time1 = (accent2 - accent1) / drum_fs; % calculate the time difference between two accents
% time2 = (accent3 - accent2) / drum_fs;
% bpm = round((1/time1 + 1/time2)*60/2) * 2; % BPM multiplied by 2 beacause the accent is on beat 1 and beat 3
% fprintf("Beats per minute: %d\n", bpm);

% Speech - detect syllables
speech_ma_best = MA_filter(x_speech, 6);
figure(plot_num);
plot_num = plot_num + 1;
start_samp = 6500;
end_samp = 11000;
plot((start_samp:end_samp), speech_ma_best(start_samp:end_samp)); % visualize input and outputs
axis([start_samp end_samp -1 1]);
xlabel("Sample number");
ylabel("Amplitude");
title("Speech - moving average filter output");


% detect number of syllables
silent = 0;
silent_counter = 0;
silent_region = [];
region_counter = 1;
lower_limit = -0.06;
upper_limit = 0.06;

for i=1:size(speech_ma_best)
    % start counting silent length
    if speech_ma_best(i) > lower_limit && speech_ma_best(i) < upper_limit
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
syllables_num = region_counter - 1;
fprintf("Number of syllables: %d\n", syllables_num);

speech_ma_best = MA_filter(x_speech, 6);
figure(plot_num);
plot((1:speech_num_x), speech_ma_best); % visualize input and outputs
axis([1 speech_num_x -1 1]);
xlabel("Sample number");
ylabel("Amplitude");
title("Speech - moving average filter output");
figure(plot_num);
for i=1:size(silent_region)
    silent_interval = silent_region(i, 2) - silent_region(i, 1);
    rectangle('Position', [silent_region(i,1),lower_limit,silent_interval,upper_limit-lower_limit]);
end
plot_num = plot_num + 1;

% detect silent region for bird
bird_ma_best = MA_filter(x_bird, 21);
silent = 0;
silent_counter = 0;
silent_region = [];
region_counter = 1;
lower_limit = -9e-3;
upper_limit = -1.4e-3;

for i=1:size(bird_ma_best)
    % start counting silent length
    if bird_ma_best(i) > lower_limit && bird_ma_best(i) < upper_limit
        silent_counter = silent_counter + 1;
    else
        if silent_counter > 2205
            silent = silent + silent_counter;
            silent_region(region_counter, :) = [i-silent_counter, i];
            region_counter = region_counter + 1;
        end
        silent_counter = 0;
    end
end
sample_time = bird_num_x / bird_fs; 
silent_time = sample_time * (silent/bird_num_x);
fprintf("Bird silent region length: %d\n", silent_time);


figure(plot_num);
plot((1:bird_num_x), bird_ma_best); 
axis([1 bird_num_x -1 1]);
xlabel("Sample number");
ylabel("Amplitude");
title("Bird - detect slient region");
figure(plot_num);
for i=1:size(silent_region)
    silent_interval = silent_region(i, 2) - silent_region(i, 1);
    rectangle('Position', [silent_region(i,1),lower_limit,silent_interval,upper_limit-lower_limit]);
end
plot_num = plot_num + 1;