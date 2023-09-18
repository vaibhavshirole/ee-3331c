%% SET MODE! 
% 
% a=1   b=2   c=3

MODE = 2;

syms s
a = 1;
b = 2;
c = 1;

if MODE == 1
    a_values = [a, a, a];
    b_values = [0.7*b, 1*b, 1.3*b];
    
    wn_results = zeros(1, 3);
    zeta_results = zeros(1, 3);
elseif MODE == 2
    a_values = [0.7*a, 1*a, 1.3*a];
    b_values = [b, b, b];
    
    wn_results = zeros(1, 3);
    zeta_results = zeros(1, 3);
elseif MODE == 3
    a_values = [a, a, a];
    b_values = [b, b, b];
    wn_values = [0.7*b, 1*b, 1.3*b];
    zeta = 0.1*a;
    
    wn_results = zeros(1, 3);
    zeta_results = zeros(1, 3);
end

for i = 1:3
    if MODE == 3
        term1 = s + zeta*wn_values(i) + wn_values(i)*sqrt(1-zeta^2)*1i;
        term2 = s + zeta*wn_values(i) - wn_values(i)*sqrt(1-zeta^2)*1i;
%     wn_result(i) = sqrt( (wn_values(i)*zeta)^2+ (wn_values(i)*sqrt(1-zeta^2))^2 );
%     zeta_result(i) = zeta;
    else
        term1 = s + a_values(i) + b_values(i)*1i;
        term2 = s + a_values(i) - b_values(i)*1i;
%     wn_values(i) = sqrt(a_values(i)^2+b_values(i)^2);
%     zeta_values(i) = a_values(i)/wn_values(i);
    end

    % Multiply the terms
    result = term1 * term2;

    % Expand the resulting polynomial
    expanded_result = expand(result);

    % Extract coefficients for x^1 and x^0
    coeffs_result = coeffs(expanded_result);
    coeff_x1 = coeffs_result(2);
    coeff_x0 = coeffs_result(1);

    % Calculate natural frequency (wn) and damping ratio (zeta)
    wn_results(i) = sqrt(coeff_x0);
    zeta_results(i) = coeff_x1 / (2 * wn_results(i));
end

% Display the results
disp('Natural Frequencies (wn):');
disp(wn_results);
disp('Damping Ratios (zeta):');
disp(zeta_results);

% Create a time vector for simulation
t = 0:0.01:10;

% Initialize a figure for plotting
figure;

% Loop through different values of wn, ζ
for i = 1:3
    zeta = zeta_results(i);
    wn = wn_results(i);
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

%% CHECK ANSWERS WITH LEC:  4-40/62
