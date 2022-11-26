% Audio preparation
function [x, x_mono, x_mono_16k, fs, num_x, new_fs, num_x_16k, plot_num] = preprocessing(infile, outfile, plot_num)
    % Preparations
    [x, fs] = audioread(infile); % read audio
    fprintf("Input - %s\n", infile);
    fprintf("  Sampling rate: %d Hz\n", fs);

    [num_x, num_channel] = size(x); % extract audio features
    fprintf("  Number of sampling point: %d\n", num_x);
    if num_channel == 1
        fprintf("  The input audio is mono\n");
        x_mono = x;
    else
        fprintf("  The input audio is stereo\n");
        x_mono = sum(x, 2)/2; 
    end

    % sound(x, fs); % play audio file

    audiowrite(outfile, x_mono, fs); % turn audio to mono and write audio to new file 

    figure(plot_num);
    plot_num = plot_num + 1;
    plot((1:num_x), x_mono); % visualize the input audio
    axis([0 num_x -1 1]);
    xlabel("Sample number");
    ylabel("Amplitude");
    plot_title = strcat(infile, " - audio input");
    title(plot_title);

    new_fs = 16000;
    x_mono_16k = resample(x_mono, new_fs, fs); % resampling
    [num_x_16k,~] = size(x_mono_16k);
    figure(plot_num);
    plot_num = plot_num + 1;
    plot((1:num_x_16k), x_mono_16k); % visualize resample audio
    axis([0 num_x_16k -1 1]);
    xlabel("Sample number");
    ylabel("Amplitude");
    plot_title = strcat(infile, " - 16kHz resampled audio");
    title(plot_title);
end