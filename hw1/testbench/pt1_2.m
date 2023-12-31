%% SET MODE! 
% 
% a=1   b=2

MODE = 2;

syms s
a = 1;
b = 2;
c = 1;

if MODE == 1
    a_values = [a, a, a];
    b_values = [0.7*b, 1*b, 1.3*b];
elseif MODE == 2
    a_values = [0.7*a, 1*a, 1.3*a];
    b_values = [b, b, b];
end

% Initialize arrays to store results
wn_values = zeros(1, 3);    % Array for natural frequencies
zeta_values = zeros(1, 3);  % Array for damping ratios

for i = 1:3
    % Define the complex conjugate terms with the fixed 'a' value and current 'b'
    term1 = s + a_values(i) + b_values(i)*1i;
    term2 = s + a_values(i) - b_values(i)*1i;

    % Multiply the terms
    result = term1 * term2;

    % Expand the resulting polynomial
    expanded_result = expand(result);

    % Extract coefficients for x^1 and x^0
    coeffs_result = coeffs(expanded_result);
    coeff_x1 = coeffs_result(2);
    coeff_x0 = coeffs_result(1);

    % Calculate natural frequency (wn) and damping ratio (zeta)
    wn_values(i) = sqrt(coeff_x0);
    zeta_values(i) = coeff_x1 / (2 * wn_values(i));

%     wn_values(i) = sqrt(a_values(i)^2+b_values(i)^2);
%     zeta_values(i) = a_values(i)/wn_values(i);
end

% Display the results
disp('Natural Frequencies (wn):');
disp(wn_values);
disp('Damping Ratios (zeta):');
disp(zeta_values);

% Create a time vector for simulation
t = 0:0.01:10;

% Initialize a figure for plotting
figure;

% Loop through different values of wn, ζ
for i = 1:3
    zeta = zeta_values(i);
    wn = wn_values(i);
    K = c * wn;

    % Calculate the closed-loop transfer function
    num = K * wn;
    den = [1, 2 * zeta * wn, wn * wn];
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
