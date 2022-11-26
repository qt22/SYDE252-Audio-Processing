% Preparations
plot_num = 1;
[drum, x_drum, x_drum_16k, drum_fs, drum_num_x, ~, ~, plot_num] = preprocessing("Drum.wav", "Drum_out_p1.wav", plot_num);
[bird, x_bird, x_bird_16k, bird_fs, bird_num_x, ~, ~, plot_num] = preprocessing("Birds.wav", "Birds_out_p1.wav", plot_num);
[speech, x_speech, x_speech_16k, speech_fs, speech_num_x, ~, ~, plot_num] = preprocessing("Speech.wav", "Speech_out_p1.wav", plot_num);

% Quantitative measures to determine window size
[drum_ma_SNR, drum_wa_SNR, drum_med_SNR, plot_num] = SNR(x_drum, "Drum", 500, 16000, 20000, 43000, 47000, plot_num);
[bird_ma_SNR, bird_wa_SNR, bird_med_SNR, plot_num] = SNR(x_bird, "Bird", 500, 85000, 87000, 89000, 91000, plot_num);
[speech_ma_SNR, speech_wa_SNR, speech_med_SNR, plot_num] = SNR(x_speech, "Speech", 500, 15000, 17000, 1, 2001, plot_num);

% Best window size
fprintf("Drum\n");
[~, drum_L_ma_best] = max(drum_ma_SNR);
fprintf("  Moving average filter, best window size: %d\n", drum_L_ma_best);
[~, drum_L_wa_best] = max(drum_wa_SNR);
fprintf("  Weighted average filter, best window size: %d\n", drum_L_wa_best);
[~, drum_L_med_best] = max(drum_med_SNR);
fprintf("Median filter, best window size: %d\n", drum_L_med_best);
figure(plot_num); plot_num = plot_num + 1;
subplot(2, 2, 1);
plot((90000:124000), x_drum(90000:124000)); % plot input
title("Drum - input");
% title("Drum - best moving average filter output");
subplot(2, 2, 2);
best_med_filter = MED_filter(x_drum, drum_L_med_best);
plot((90000:124000), best_med_filter(90000:124000)); % median filter plot
title("Drum - Median Filter");
subplot(2, 2, 3);
best_wa_filter = WA_filter(x_drum, drum_L_wa_best);
plot((90000:124000), best_wa_filter(90000:124000)); % weighted moving average filter plot
title("Drum - Weighted Moving Average Filter");
subplot(2, 2, 4);
best_ma_filter = MA_filter(x_drum, drum_L_ma_best);
plot((90000:124000), best_ma_filter(90000:124000)); % moving average filter plot
title("Drum - Moving Average Filter");
% Drum best filter:
drum_best = best_ma_filter;


fprintf("Bird\n");
[~, bird_L_ma_best] = max(bird_ma_SNR); % SNR is not a good metric for bird
bird_L_ma_best = 21; % determined by subjective measurements
fprintf("  Moving average filter, best window size: %d\n", bird_L_ma_best);
[~, bird_L_wa_best] = max(bird_wa_SNR); % SNR is not a good metric for bird
bird_L_wa_best = 9; % determined by subjective measurements
fprintf("  Weighted average filter, best window size: %d\n", bird_L_wa_best);
[~, bird_L_med_best] = max(bird_med_SNR); % SNR is not a good metric for bird
bird_L_med_best = 2; % determined by subjective measurements
fprintf("Median filter, best window size: %d\n", bird_L_med_best);
% Bird best filter:
bird_best = MA_filter(x_bird, bird_L_ma_best);

fprintf("Speech\n");
[~, speech_L_ma_best] = max(speech_ma_SNR);
fprintf("  Moving average filter, best window size: %d\n", speech_L_ma_best);
[~, speech_L_wa_best] = max(speech_wa_SNR);
fprintf("  Weighted average filter, best window size: %d\n", speech_L_wa_best);
[~, speech_L_med_best] = max(speech_med_SNR);
fprintf("Median filter, best window size: %d\n", speech_L_med_best);
% subplots
figure(plot_num); plot_num = plot_num + 1;
subplot(2, 2, 1);
plot((12000:22000), x_speech(12000:22000)); % plot input
title("Speech - input");
subplot(2, 2, 2);
best_med_filter = MED_filter(x_speech, speech_L_med_best);
plot((12000:22000), best_med_filter(12000:22000)); % median filter plot
title("Speech - Median Filter");
subplot(2, 2, 3);
best_wa_filter = WA_filter(x_speech, speech_L_wa_best);
plot((12000:22000), best_wa_filter(12000:22000)); % weighted moving average filter plot
title("Speech - Weighted Moving Average Filter");
subplot(2, 2, 4);
best_ma_filter = MA_filter(x_speech, speech_L_ma_best);
plot((12000:22000), best_ma_filter(12000:22000)); % moving average filter plot
title("Speech - Moving Average Filter");
% Speech best filter:
speech_best = best_wa_filter;

% Filtering
drum_ma_slight_distort = MA_filter(x_drum, 13);
drum_ma_very_distort = MA_filter(x_drum, 890);
figure(plot_num);
plot_num = plot_num + 1;
plot((1:drum_num_x), x_drum); % visualize input
axis([0 drum_num_x -1 1]);
xlabel("Sample number");
ylabel("Amplitude");
title("Drum - moving average filter output");
hold on;
plot((1:drum_num_x), drum_ma_slight_distort); % overlay slightly distorted window size plot
plot((1:drum_num_x), drum_ma_very_distort); % overlay very distorted window size plot
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
plot((1:bird_num_x), bird_ma_slight_distort); % overlay slightly distorted window size plot
plot((1:bird_num_x), bird_ma_very_distort); % overlay very distorted window size plot
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
plot((1:speech_num_x), speech_ma_slight_distort); % overlay slightly distorted window size plot
plot((1:speech_num_x), speech_ma_very_distort); % overlay very distorted window size plot
legend('Input', 'Window Size: 5 (Slightly distorted)', 'Window Size: 1500 (Very distorted)');
hold off;

% Weighted average filter
drum_wa_slight_distort = WA_filter(x_drum, 5);
drum_wa_very_distort = WA_filter(x_drum, 1477);
figure(plot_num);
plot_num = plot_num + 1;
plot((1:drum_num_x), x_drum); % visualize input and outputs
axis([0 drum_num_x -1 1]);
xlabel("Sample number");
ylabel("Amplitude");
title("Drum - weighted average filter output");
hold on;
plot((1:drum_num_x), drum_wa_slight_distort); % overlay slightly distorted window size plot
plot((1:drum_num_x), drum_wa_very_distort); % overlay very distorted window size plot
legend('Input', 'Window Size: 5 (Slightly distorted)', 'Window Size: 1477 (Very distorted)');
hold off;

bird_wa_slight_distort = WA_filter(x_bird, 5);
bird_wa_very_distort = WA_filter(x_bird, 217);
figure(plot_num);
plot_num = plot_num + 1;
plot((1:bird_num_x), x_bird); % visualize input and outputs
axis([0 bird_num_x -1 1]);
xlabel("Sample number");
ylabel("Amplitude");
title("Bird - weighted average filter output");
hold on;
plot((1:bird_num_x), bird_wa_slight_distort); % overlay slightly distorted window size plot
plot((1:bird_num_x), bird_wa_very_distort); % overlay very distorted window size plot
legend('Input', 'Window Size: 5 (Slightly distorted)', 'Window Size: 217 (Very distorted)');
hold off;

speech_wa_slight_distort = WA_filter(x_speech, 7);
speech_wa_very_distort = WA_filter(x_speech, 1310);
figure(plot_num);
plot_num = plot_num + 1;
plot((1:speech_num_x), x_speech); % visualize input and outputs
axis([0 speech_num_x -1 1]);
xlabel("Sample number");
ylabel("Amplitude");
title("Speech - weighted average filter output");
hold on;
plot((1:speech_num_x), speech_wa_slight_distort); % overlay slightly distorted window size plot
plot((1:speech_num_x), speech_wa_very_distort); % overlay very distorted window size plot
legend('Input', 'Window Size: 7 (Slightly distorted)', 'Window Size: 1310 (Very distorted)');
hold off;


% Median filter
drum_med_slight_distort = MED_filter(x_drum, 4);
drum_med_very_distort = MED_filter(x_drum, 398);
figure(plot_num);
plot_num = plot_num + 1;
plot((1:drum_num_x), x_drum); % visualize input and outputs
axis([0 drum_num_x -1 1]);
xlabel("Sample number");
ylabel("Amplitude");
title("Drum - median filter output");
hold on;
plot((1:drum_num_x), drum_med_slight_distort); % overlay slightly distorted window size plot
plot((1:drum_num_x), drum_med_very_distort); % overlay very distorted window size plot
legend('Input', 'Window Size: 4 (Slightly distorted)', 'Window Size: 398 (Very distorted)');
hold off;

bird_med_slight_distort = MA_filter(x_bird, 2);
bird_med_very_distort = MA_filter(x_bird, 37);
figure(plot_num);
plot_num = plot_num + 1;
plot((1:bird_num_x), x_bird); % visualize input and outputs
axis([0 bird_num_x -1 1]);
xlabel("Sample number");
ylabel("Amplitude");
title("Bird - median filter output");
hold on;
plot((1:bird_num_x), bird_med_slight_distort); % overlay slightly distorted window size plot
plot((1:bird_num_x), bird_med_very_distort); % overlay very distorted window size plot
legend('Input', 'Window Size: 2 (Slightly distorted)', 'Window Size: 37 (Very distorted)');
hold off;

speech_med_slight_distort = MA_filter(x_speech, 3);
speech_med_very_distort = MA_filter(x_speech, 37);
figure(plot_num);
plot_num = plot_num + 1;
plot((1:speech_num_x), x_speech); % visualize input and outputs
axis([0 speech_num_x -1 1]);
xlabel("Sample number");
ylabel("Amplitude");
title("Speech - median filter output");
hold on;
plot((1:speech_num_x), speech_med_slight_distort); % overlay slightly distorted window size plot
plot((1:speech_num_x), speech_med_very_distort); % overlay very distorted window size plot
legend('Input', 'Window Size: 3 (Slightly distorted)', 'Window Size: 37 (Very distorted)');
hold off;

% Analyzing signals
% Drum - detect bpm
interval = 4000;
acc1_start = 16000;
acc1_end = acc1_start + interval;
acc2_start = 46000;
acc2_end = acc2_start + interval;
acc3_start = 78000;
acc3_end = acc3_start + interval;
[~, accent1] = max(drum_best(acc1_start:acc1_end)); % find the peak within the time of interest
accent1 = accent1 + acc1_start - 1;
[~, accent2] = max(drum_best(acc2_start:acc2_end));
accent2 = accent2 + acc2_start - 1;
[~, accent3] = max(drum_best(acc3_start:acc3_end));
accent3 = accent3 + acc3_start - 1;
time1 = (accent2 - accent1) / drum_fs; % calculate the time difference between two accents
time2 = (accent3 - accent2) / drum_fs;
bpm = round((1/time1 + 1/time2)*60/2*2); % BPM multiplied by 2 beacause the accent is on beat 1 and beat 3
fprintf("Beats per minute: %12.12f\n", bpm);

figure(plot_num);
plot((1:drum_num_x), drum_best); % plot drum moving average best filter
axis([1 drum_num_x -1 1]);
xlabel("Sample number");
ylabel("Amplitude");
title("Drum - detect BPM");
figure(plot_num);
rectangle('Position', [acc1_start,-0.2,interval,0.4]);
rectangle('Position', [acc2_start,-0.2,interval,0.4]);
rectangle('Position', [acc3_start,-0.2,interval,0.4]);
plot_num = plot_num + 1;


% Speech - detect syllables
pause_counter = 0;
pause_region = [];
region_counter = 1;
lower_limit = -0.06;
upper_limit = 0.06;
min_pause_length = 1000;

for i=1:size(speech_best)
    % iterate through the whole sample to locate quiet samples
    if speech_best(i) > lower_limit && speech_best(i) < upper_limit 
        pause_counter = pause_counter + 1;
    else
        % if the pause samples last longer than the minimum defined range
        % (0.0625 seconds)
        if pause_counter > min_pause_length
            % add the start and end index to a Lx2 matrix to record all pause regions
            pause_region(region_counter, :) = [i-pause_counter, i];
            region_counter = region_counter + 1;
        end
        pause_counter = 0;
    end

    % check if last part is a pause
    if i == speech_num_x
        if pause_counter > min_pause_length
            % add last part of the signal to the pause region matrix
            pause_region(region_counter, :) = [i-pause_counter, i];
        end
    end
end
syllables_num = region_counter - 1; % syllable count is one less than the pause regions
fprintf("Number of syllables: %d\n", syllables_num);

figure(plot_num);
plot((1:speech_num_x), speech_best); % plot best speech filter
axis([1 speech_num_x -1 1]);
xlabel("Sample number");
ylabel("Amplitude");
title("Speech - weighted moving average filter output");
figure(plot_num);
for i=1:size(pause_region)
    silent_interval = pause_region(i, 2) - pause_region(i, 1);
    % box the pause regions
    rectangle('Position', [pause_region(i,1),lower_limit,silent_interval,upper_limit-lower_limit]);
end
plot_num = plot_num + 1;



% Bird - detect silent regions 
silent = 0;
silent_counter = 0;
silent_region = [];
region_counter = 1;
lower_limit = -9e-3;
upper_limit = -1.4e-3;
min_silent_length = 2205;

for i=1:size(bird_best)
     % iterate through the whole sample to locate quiet samples
    if bird_best(i) > lower_limit && bird_best(i) < upper_limit
        silent_counter = silent_counter + 1;
    else
        % if the silent samples last longer than the minimum defined range
        % (0.2 seconds)
        if silent_counter > min_silent_length
            % add the start and end index to a Lx2 matrix to record all pause regions
            silent = silent + silent_counter;
            silent_region(region_counter, :) = [i-silent_counter, i];
            region_counter = region_counter + 1;
        end
        silent_counter = 0;
    end

     % check if last part is a pause
    if i == bird_num_x
        if silent_counter > min_silent_length
            % add last part of the signal to the pause region matrix
            silent_region(region_counter, :) = [i-silent_counter, i];
        end
    end
end
sample_time = bird_num_x / bird_fs; % audio file length 
silent_time = sample_time * (silent/bird_num_x); % audio file length * silent ratio = silent time
fprintf("Bird silent region length: %d\n", silent_time);


figure(plot_num);
plot((1:bird_num_x), bird_best); % plot best bird filter
axis([1 bird_num_x -1 1]);
xlabel("Sample number");
ylabel("Amplitude");
title("Bird - Silent Region");
figure(plot_num);
for i=1:size(silent_region)
    silent_interval = silent_region(i, 2) - silent_region(i, 1);
    % box the silent regions
    rectangle('Position', [silent_region(i,1),lower_limit,silent_interval,upper_limit-lower_limit]);
end
plot_num = plot_num + 1;