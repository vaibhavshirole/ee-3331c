% Define system parameters
a = 1;
b = 2;
c = 1;
K = c;

% Define different values for wn
%wn_values = [74^(1/2)/5, 5^(1/2), 194^(1/2)/5];    %partA
wn_values = [449^(1/2)/10, 5^(1/2), 569^(1/2)/10];  %partB

% Define different values of damping ratio ζ
%zeta_values = [(5*74^(1/2))/74, 5^(1/2)/5, (5*194^(1/2))/194];
zeta_values = [(7*449^(1/2))/449, 5^(1/2)/5, (13*569^(1/2))/569];

% Create a time vector for simulation
t = 0:0.01:10;

% Initialize a figure for plotting
figure;

% Loop through different values of wn, ζ
for i = 1:length(zeta_values)
    zeta = zeta_values(i);
    wn = wn_values(i);

    % Calculate the poles
    s1 = -a + 1j*b;
    s2 = -a - 1j*b;

    % Calculate the closed-loop transfer function
    num = K * wn;
    den = [1, 2 * zeta * wn, wn^2];
    sys = tf(num, den);

    % Simulate the step response
    y = step(sys, t);

    % Plot the step response
    plot(t, y, 'LineWidth', 1.5, 'DisplayName', ['ζ = ', num2str(zeta)]);
    hold on;
end

% Add labels and legend
xlabel('Time (s)');
ylabel('Amplitude');
title('Step Response of Second-Order System with Varying ζ');
legend('Location', 'best');
grid on;

% Display the plot
hold off;
