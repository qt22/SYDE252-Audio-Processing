% Median filter
function y = MED_filter(x, L)
    [sample_num, ~] = size(x); 
    y = zeros(size(x)); % initialize output with 0s
    x = cat(1,zeros([L-1, 1]),x); % add zeros in front of x to compensate for the first element
    for n = 1:sample_num
        y(n) = median(x(n:n+L-1));
    end
end