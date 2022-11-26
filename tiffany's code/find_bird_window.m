[bird, fs] = audioread("Birds.wav"); 
[bird_num_x, ~] = size(bird);

plot_num = 1;

for i=21:2:21
    disp("window size: " + i);
    sound(MA_filter(bird, i)*i, fs); % sound(MA_filter(talk, 500)*500, fs)
    pause(10.5);
    
    figure(plot_num);
    plot_num = plot_num + 1;
    plot((1:bird_num_x), bird); % visualize input
    axis([0 bird_num_x -1 1]);
    xlabel("Sample number");
    ylabel("Amplitude");
    title("Bird - moving average filter output");
    
    hold on;
    plot((1:bird_num_x), MA_filter(bird, i)*i); 
    legend('Input', 'Moving Average Filter Output');
    hold off;
end