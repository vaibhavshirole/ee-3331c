syms s
a = 1;
b = 2;

% Values of 'a' for three iterations
a_values = [0.7*a, 1*a, 1.3*a];

% Initialize arrays to store results
wn_values = zeros(1, 3);    % Array for natural frequencies
zeta_values = zeros(1, 3);  % Array for damping ratios

for i = 1:3
    % Define the complex conjugate terms with the current 'a' value and fixed 'b'
    term1 = s + a_values(i) + b*1i;
    term2 = s + a_values(i) - b*1i;

    % Multiply the terms
    result = term1 * term2;

    % Expand the resulting polynomial
    expanded_result = expand(result);

    % Extract coefficients for x^1 and x^0
    coeffs_result = coeffs(expanded_result);

    % Extract the coefficient for x^1 (s, 1)
    coeff_x1 = coeffs_result(2);

    % Extract the coefficient for x^0 (s, 0)
    coeff_x0 = coeffs_result(1);

    % Calculate natural frequency (wn) and damping ratio (zeta)
    wn = sqrt(coeff_x0);
    zeta = coeff_x1 / (2 * wn);

    % Store results in arrays
    wn_values(i) = wn;
    zeta_values(i) = zeta;
end

% Display the results
disp('Natural Frequencies (wn):');
disp(wn_values);
disp('Damping Ratios (zeta):');
disp(zeta_values);