% Weighted average filter
function y = WA_filter(x, L)
    b = gausswin(L); % generate guase
    b_norm = b / sum(b);
    y = filter(b_norm, 1, x);
    % wvtool(gausswin(L));
end