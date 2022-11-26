[speech, fs1] = audioread("Speech.wav"); 
[drum, ~] = audioread("Drum.wav");
[drum_num_x, ~] = size(drum);
[bird, ~] = audioread("Birds.wav");

plot_num = 1;

drum_ma_slight_distort = MA_filter(drum, 13);
% drum_ma_very_distort = MA_filter(drum, 500);
figure(plot_num);
plot_num = plot_num + 1;
plot((1:drum_num_x), drum); % visualize input
axis([0 drum_num_x -1 1]);
xlabel("Sample number");
ylabel("Amplitude");
title("Drum - moving average filter output");
hold on;
plot((1:drum_num_x), drum_ma_slight_distort); % slight distort 
% plot((1:drum_num_x), drum_ma_very_distort);
legend('Input', 'Slightly distorted filter output');
hold off;