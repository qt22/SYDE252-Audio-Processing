% Preparations
drum_infile = "Drum.wav";
drum_outfile = "Drum_out_p1.wav";
[x, fs] = audioread(drum_infile); % read audio
fprintf("Sampling rate: %d Hz\n", fs);

[num_x, num_channel] = size(x); % extract audio features
fprintf("Number of sampling point: %d\n", num_x);
if num_channel == 1
    fprintf("The input audio is mono\n");
    x_mono = x;
else
    fprintf("The input audio is stereo\n");
    x_mono = sum(x, 2)/2; 
end

plot_num = 1;

% sound(x_mono, fs); % play audio file

audiowrite(drum_outfile, x_mono, fs); % turn audio to mono and write audio to new file 

figure(plot_num);
plot_num = plot_num + 1;
plot((1:num_x), x_mono); % visualize the input audio
axis([0 num_x -1 1]);
xlabel("Sample number");
ylabel("Amplitude");
title("Drum audio input");

% new_fs = 16000;
% x_mono_16k = resample(x_mono, new_fs, fs); % resampling
% [num_x_16k,~] = size(x_mono_16k);
% figure(plot_num);
% plot_num = plot_num + 1;
% plot((1:num_x_16k), x_mono_16k); % visualize resample audio
% axis([0 num_x_16k -1 1]);
% xlabel("Sample number");
% ylabel("Amplitude");
% title("16kHz resampled drum audio");

% Quantitative measures to determine window size
% Signal - sample number: 6000 - 8000
% Noise - sample number: 16400 - 17400
signal_start = 16000;
signal_end = 20000;
noise_start = 43000;
noise_end = 47000;
x_SNR = 20*log10(sum(x_mono(signal_start:signal_end).*x_mono(signal_start:signal_end))/sum(x_mono(noise_start:noise_end).*x_mono(noise_start:noise_end)));
fprintf("Input SNR: %12.12f\n", x_SNR);
L_MAX = 500;
ma_SNR = zeros([L_MAX, 1]);
wa_SNR = zeros([L_MAX, 1]);
med_SNR = zeros([L_MAX, 1]);
for L = 1:L_MAX
    y_ma = MA_filter(x_mono, L)*L;
    y_wa = WA_filter(x_mono, L)*L;
    y_med = MED_filter(x_mono, L)*L;

    ma_SNR(L) = 20*log10(sum(y_ma(signal_start:signal_end).*y_ma(signal_start:signal_end))/sum(y_ma(noise_start:noise_end).*y_ma(noise_start:noise_end)));
    wa_SNR(L) = 20*log10(sum(y_wa(signal_start:signal_end).*y_wa(signal_start:signal_end))/sum(y_wa(noise_start:noise_end).*y_wa(noise_start:noise_end)));
    med_SNR(L) = 20*log10(sum(y_med(signal_start:signal_end).*y_med(signal_start:signal_end))/sum(y_med(noise_start:noise_end).*y_med(noise_start:noise_end)));
end

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

% Filtering
% Moving average filter
[~, L_ma_best] = max(ma_SNR);
fprintf("Moving average filter, best window size: %d\n", L_ma_best);
y_ma_best = MA_filter(x_mono, L_ma_best)*20; % introduce a gain of 20 to compensate the attenuation
y_ma_best = max(y_ma_best, -1);
y_ma_best = min(y_ma_best, 1);
ma_noise = abs(y_ma_best-x_mono);
figure(plot_num);
plot_num = plot_num + 1;
plot((1:num_x), y_ma_best); % visualize input and outputs
axis([0 num_x -1 1]);
xlabel("Sample number");
ylabel("Amplitude");
title("Moving average filter output");

% Weighted average filter
[~, L_wa_best] = max(wa_SNR); % output is noisy, require subjective analysis
L_wa_best = 30; % value obtained by subjective analysis
fprintf("Weighted average filter, best window size: %d\n", L_wa_best);
y_wa_best = WA_filter(x_mono, L_wa_best)*20;
y_wa_best = max(y_wa_best, -1);
y_wa_best = min(y_wa_best, 1);
figure(plot_num);
plot_num = plot_num + 1;
plot((1:num_x), y_wa_best); % visualize input and outputs
axis([0 num_x -1 1]);
xlabel("Sample number");
ylabel("Amplitude");
title("Weighted average filter output");

% Median filter
[~, L_med_best] = max(med_SNR);
fprintf("Median filter, best window size: %d\n", L_med_best);
y_med_best = MED_filter(x_mono, L_med_best)*10;
y_med_best = max(y_med_best, -1);
y_med_best = min(y_med_best, 1);
figure(plot_num);
plot_num = plot_num + 1;
plot((1:num_x), y_med_best);
axis([0 num_x -1 1]);
xlabel("Sample number");
ylabel("Amplitude");
title("Median filter output");

% Analyzing signals
acc1_start = 16000;
acc1_end = 20000;
acc2_start = 46000;
acc2_end = 50000;
acc3_start = 78000;
acc3_end = 82000;
[~, accent1] = max(y_ma_best(acc1_start:acc1_end)); % find the peak within the time of interest
accent1 = accent1 + acc1_start - 1;
[~, accent2] = max(y_ma_best(acc2_start:acc2_end));
accent2 = accent2 + acc2_start - 1;
[~, accent3] = max(y_ma_best(acc3_start:acc3_end));
accent3 = accent3 + acc3_start - 1;
time1 = (accent2 - accent1) / fs; % calculate the time difference between two accents
time2 = (accent3 - accent2) / fs;
bpm = round((1/time1 + 1/time2)*60/2) * 2; % BPM multiplied by 2 beacause the accent is on beat 1 and beat 3
fprintf("Beats per minute: %d\n", bpm);

% % Evaluate noise after filtering
% % FFT of input - selected region
% fft_fs = 8000;
% start_samp = 28500;
% end_samp = 29500; % num_x;
% length = end_samp - start_samp;
% freq = fft_fs*(0:floor(length/2))/length; % Get rid of mirrored output
% x_fft = fft(x_mono(start_samp:end_samp));
% x_fft = mag2db(abs(x_fft));
% x_fft = x_fft((1:(floor(length/2))+1));
% figure(plot_num);
% plot_num = plot_num + 1;
% plot(freq, x_fft);
% xlabel("Frequency (Hz)");
% ylabel("Magnitude (dB)");
% title("FFT of a portion of the audio input");
% 
% % FFT of input - entire clip
% length = num_x;
% freq = fft_fs*(0:floor(length/2))/length; % Get rid of mirrored output
% x_fft = fft(x_mono);
% x_fft = mag2db(abs(x_fft));
% x_fft = x_fft((1:(floor(length/2))+1));
% figure(plot_num);
% plot_num = plot_num + 1;
% plot(freq, x_fft);
% xlabel("Frequency (Hz)");
% ylabel("Magnitude (dB)");
% title("FFT of audio input");
% 
% % FFT of the output of moving average filter
% freq = fft_fs*(0:floor(length/2))/length; % Get rid of mirrored output
% ma_fft = fft(MA_filter(x_mono,12)*12);
% ma_fft = mag2db(abs(ma_fft));
% ma_fft = ma_fft((1:(floor(length/2))+1));
% figure(plot_num);
% plot_num = plot_num + 1;
% plot(freq, ma_fft);
% xlabel("Frequency (Hz)");
% ylabel("Magnitude (dB)");
% title("FFT of moving average output");
% figure(plot_num);
% plot_num = plot_num + 1;
% plot(freq, (ma_fft - x_fft));
% xlabel("Frequency (Hz)");
% ylabel("Magnitude (dB)");
% title("Frequency response of moving average filter");
% 
% % FFT of the output of weighted average filter
% freq = fft_fs*(0:floor(length/2))/length; % Get rid of mirrored output
% wa_fft = fft(y_wa_best);
% wa_fft = mag2db(abs(wa_fft));
% wa_fft = wa_fft((1:(floor(length/2))+1));
% figure(plot_num);
% plot_num = plot_num + 1;
% plot(freq, wa_fft);
% xlabel("Frequency (Hz)");
% ylabel("Magnitude (dB)");
% title("FFT of weighted average output");
% figure(plot_num);
% plot_num = plot_num + 1;
% plot(freq, (wa_fft - x_fft));
% xlabel("Frequency (Hz)");
% ylabel("Magnitude (dB)");
% title("Frequency response of weighted average filter");
% 
% % FFT of the output of median filter
% freq = fft_fs*(0:floor(length/2))/length; % Get rid of mirrored output
% med_fft = fft(y_med_best);
% med_fft = mag2db(abs(med_fft));
% med_fft = med_fft((1:(floor(length/2))+1));
% figure(plot_num);
% plot_num = plot_num + 1;
% plot(freq, med_fft);
% xlabel("Frequency (Hz)");
% ylabel("Magnitude (dB)");
% title("FFT of median output");
% figure(plot_num);
% plot_num = plot_num + 1;
% plot(freq, (med_fft - x_fft));
% xlabel("Frequency (Hz)");
% ylabel("Magnitude (dB)");
% title("Frequency response of median filter");