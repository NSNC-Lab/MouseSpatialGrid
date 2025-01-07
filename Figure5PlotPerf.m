% Example data
categories = {'SPIKE', 'ISI', 'RISPIKE'}; % Labels for categories
performance = [mean(perf.SPIKE),mean(perf.ISI), mean(perf.RISPIKE)]; % Accuracy values in percentage

% Create a bar graph
figure;
bar(performance, 'FaceColor', [0.2 0.6 0.8]); % Create bar graph with custom color

% Add labels and title
set(gca, 'XTickLabel', categories); % Set X-axis labels
xlabel('Models'); % Label for X-axis
ylabel('Performance'); % Label for Y-axis
title('Model Accuracy Comparison'); % Title

% Add grid lines for better readability
grid on;

% Annotate bars with their values
for i = 1:length(accuracy)
    text(i, accuracy(i) + 2, sprintf('%d%%', accuracy(i)), ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom');
end

