% Transfer function G(s)
numerator = conv([10, 1], [1, 1]); % (10s + 1)(s + 1)
denominator = conv([100, 1], [0.1, 1]); % (100s + 1)(0.1s + 1)

G = tf(numerator, denominator);

% Bode plot to visualize the frequency response
figure;
bode(G);
grid on;
title('Bode Plot of G(s)');

