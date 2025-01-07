% Parameters
t = 0:0.01:10; % Time vector
tau_I = 0.06;     % Time constant for w_I
tau_A = 0.1;     % Time constant for w_A

% Input signal: Step function
s = double(t > 2 & t < 8);

% Kernel definitions (exponential decay)
w_I = exp(-t/tau_I) .* (t >= 0);
w_A = exp(-t/tau_A) .* (t >= 0);

% Normalize kernels
w_I = w_I / sum(w_I);
w_A = w_A / sum(w_A);

% Step 1: Intensity Integration
r_I = conv(s, w_I, 'same'); % Convolve with w_I

% Step 2: Adaptive Gain
adaptive_term = conv(r_I, w_A, 'same'); % Convolve with w_A
g = 1 ./ (1 + adaptive_term); % Adaptive gain

% Step 3: Intensity Gain Control
r_IA = g .* r_I;

% Plot Results
figure;
subplot(4, 1, 1);
plot(t, s, 'k', 'LineWidth', 1.5);
title('Input Signal s(t)');
xlabel('Time');
ylabel('Amplitude');

subplot(4, 1, 2);
plot(t, r_I, 'b', 'LineWidth', 1.5);
title('Intensity Integration r_I(t)');
xlabel('Time');
ylabel('Amplitude');

subplot(4, 1, 3);
plot(t, g, 'r', 'LineWidth', 1.5);
title('Adaptive Gain g(t)');
xlabel('Time');
ylabel('Gain');

subplot(4, 1, 4);
plot(t, r_IA, 'g', 'LineWidth', 1.5);
title('Intensity Gain Control r_{IA}(t)');
xlabel('Time');
ylabel('Amplitude');