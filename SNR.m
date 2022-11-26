% SNR - quantitative measures to determine window size
function [ma_SNR, wa_SNR, med_SNR, plot_num] = SNR(x, clip, L_max, signal_start, signal_end, noise_start, noise_end, plot_num)
    x_SNR = 20*log10(sum(x(signal_start:signal_end).*x(signal_start:signal_end))/sum(x(noise_start:noise_end).*x(noise_start:noise_end)));
    fprintf("%s - Input SNR: %12.12f\n", clip, x_SNR);
    ma_SNR = zeros([L_max, 1]);
    wa_SNR = zeros([L_max, 1]);
    med_SNR = zeros([L_max, 1]);
    for L = 1:L_max
        y_ma = MA_filter(x, L)*L;
        y_wa = WA_filter(x, L)*L;
        y_med = MED_filter(x, L)*L;

        ma_SNR(L) = 20*log10(sum(y_ma(signal_start:signal_end).*y_ma(signal_start:signal_end))/sum(y_ma(noise_start:noise_end).*y_ma(noise_start:noise_end)));
        wa_SNR(L) = 20*log10(sum(y_wa(signal_start:signal_end).*y_wa(signal_start:signal_end))/sum(y_wa(noise_start:noise_end).*y_wa(noise_start:noise_end)));
        med_SNR(L) = 20*log10(sum(y_med(signal_start:signal_end).*y_med(signal_start:signal_end))/sum(y_med(noise_start:noise_end).*y_med(noise_start:noise_end)));
    end

    figure(plot_num);
    plot_num = plot_num + 1;
    plot(1:L_max, ma_SNR);
    hold on;
    plot(1:L_max, wa_SNR);
    plot(1:L_max, med_SNR);
    legend('Moving average filter', 'Weighted average filter', 'Median filter');
    xlabel("Window size");
    ylabel("SNR (dB)");
    plot_title = strcat(clip, " - SNR ratio");
    title(plot_title);
    hold off;
end