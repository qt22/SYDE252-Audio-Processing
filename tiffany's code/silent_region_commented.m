%%% part 3: signal analysis
% detect silent region

% silent: total sample counter for all determined silent regions
silent = 0;

% temp counter to determine the sample length of each potential silent
% region
silent_counter = 0;


silent_region = [];

% counts total number of silent regions 
region_counter = 1;

% define acceptable range with min and max amplitude of silent regions
lower_limit = -5.85e-3;
upper_limit = -3.8e-3;

% iterate from the first sample to the last sample of filtered output
for i=1:size(y_ma)

    % if sample is within acceptable "silent amplitude"
    if y_ma(i) > lower_limit && y_ma(i) < upper_limit
        % add 1 to the sample counter 
        silent_counter = silent_counter + 1;

    % if sample outside of acceptable "silent amplitude" (birds chirping)
    else
        % check if silent region was greater than the min sample len
        if silent_counter > 3000
            
            % add samples to the total counter for silent region
            silent = silent + silent_counter;

            % record the start and end indexs of the silent region in an L x 2 matrix
            silent_region(region_counter, :) = [i-silent_counter, i];

            % keep track of number of silent regions identified 
            region_counter = region_counter + 1;
        end
        % reset temp counter
        silent_counter = 0;
    end
end
% number of samples in input signal
[sample_size, ~] = size(bird);
% total time of signal
sample_time = sample_size / fs;
% find the total seconds of silent time using ratio of silent samples over
% total samples
silent_time = sample_time * (silent/sample_size);
disp("silent time: " + silent_time);

% plot boxes around the silent regions
figure(41);
for i=1:size(silent_region)
    silent_interval = silent_region(i, 2) - silent_region(i, 1);
    rectangle('Position', [silent_region(i,1),lower_limit,silent_interval,upper_limit-lower_limit]);
end
