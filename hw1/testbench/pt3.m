syms s wn
a = 1;
b = 2;
c = 1;
K = c;

% Values for 'wn' in three iterations
wn_values = [0.7*b, 1*b, 1.3*b];
zeta = 0.1*a;

% Initialize arrays to store results
wn_result = zeros(1, 3);    % Array for natural frequencies
zeta_result = zeros(1, 3);  % Array for damping ratios

% Loop through three different values of 'wn'
for i = 1:3
    % Define the complex conjugate terms with the current 'a' value and fixed 'b'
    term1 = s + zeta*wn_values(i) + wn_values(i)*sqrt(1-zeta^2)*1i;
    term2 = s + zeta*wn_values(i) - wn_values(i)*sqrt(1-zeta^2)*1i;

    % Multiply the terms
    result = term1 * term2;

    % Expand the resulting polynomial
    expanded_result = expand(result);

    % Extract coefficients for x^1 and x^0
    coeffs_result = coeffs(expanded_result);
    coeff_x1 = coeffs_result(2);
    coeff_x0 = coeffs_result(1);

    % Calculate natural frequency (wn) and damping ratio (zeta)
    wn_result(i) = sqrt(coeff_x0);
    zeta_result(i) = coeff_x1 / (2 * wn_result(i));

%     wn_result(i) = sqrt( (wn_values(i)*zeta)^2+ (wn_values(i)*sqrt(1-zeta^2))^2 );
%     zeta_result(i) = zeta;
end

% Display the results
disp('Natural Frequencies (wn):');
disp(wn_result);
disp('Damping Ratios (zeta):');
disp(zeta_result);

% Create a time vector for simulation (low wn, increased time)
t = 0:0.01:25;

% Initialize a figure for plotting
figure;

% Loop through different values of wn, ζ
for i = 1:3
    zeta = zeta_result(i);
    wn = wn_result(i);
    K = c * wn;

    % Calculate the closed-loop transfer function
    num = K * wn;
    den = [1, 2 * zeta * wn, wn * wn];
    sys = tf(num, den);

    % Simulate the step response
    y = step(sys, t);

    % Plot the step response
    plot(t, y, 'LineWidth', 1.5, 'DisplayName', ['wn = ', num2str(wn)]);
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

