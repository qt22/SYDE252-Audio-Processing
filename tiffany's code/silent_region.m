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
