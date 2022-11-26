% MA
function y = MA_filter(x, L)
    b = (1/L)*ones(1,L);
    y = filter(b, 1, x);
%     y = x - y; % HPF
end