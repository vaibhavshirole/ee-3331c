%% A //----------//----------//----------//----------//----------//

a = 1;
b = 2;
b_values = [0.05*b, 1*b, 50*b];

figure;
for i = 1:length(b_values)
    b = b_values(i);
    G = tf(a, [1, b])
    
    % Calculate Bode plot
    [m, p, w] = bode(G);
    
    % Convert angular freq to freq
    f = w/(2*pi);
    
    % Plot magnitude response
    subplot(2, 1, 1);
    semilogx(f, 20*log10(squeeze(m)), 'DisplayName', ['b = ' num2str(b)], 'LineWidth', 2);
    hold on;
    
    % Plot phase response
    subplot(2, 1, 2);
    semilogx(f, squeeze(p), 'DisplayName', ['b = ' num2str(b)], 'LineWidth', 2);
    hold on;
end

% Set labels and legends and show plot
subplot(2, 1, 1);
grid on;
title('Bode Plot - Magnitude Response');
xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)');
legend('Location', 'Best');

subplot(2, 1, 2);
grid on;
title('Bode Plot - Phase Response');
xlabel('Frequency (Hz)'); ylabel('Phase (degrees)');
legend('Location', 'Best');

sgtitle('Bode Plots for Different Values of b');

%% B //----------//----------//----------//----------//----------//

a = 1;
b = 2;
b_values = [0.05*b, 1*b, 50*b];

figure;
for i = 1:length(b_values)
    b = b_values(i);
    G = tf(a, [1, b])
    
    % Calculate the Bode plot
    [m, p] = bode(G);
    
    % Make into 1D array
    m = m(1, :);
    p = p(1, :);
    
    % Convert to radians 
    polarplot((p*pi/180), m, 'DisplayName', ['b = ' num2str(b)], 'LineWidth', 2);
    hold on;
end
legend('Location', 'Best');

%% C //----------//----------//----------//----------//----------//

b = 2;
c_values = [0.1*b, 10*b];

% Magnitude and phase plots
figure;
for i = 1:length(c_values)
    c = c_values(i);
    G = tf([b, 1], [c, 1, 0, 0])
    
    % Calculate the Bode plot
    [m, p, w] = bode(G);
    
    % Convert angular freq to freq
    f = w/(2*pi);
    
    % Plot magnitude response
    subplot(2, 1, 1);
    semilogx(f, 20*log10(squeeze(m)), 'DisplayName', ['c = ' num2str(c)], 'LineWidth', 2);
    hold on;
    
    % Plot phase response
    subplot(2, 1, 2);
    semilogx(f, squeeze(p), 'DisplayName', ['c = ' num2str(c)], 'LineWidth', 2);
    hold on;
end

% Set labels and legends and show plot
subplot(2, 1, 1);
grid on;
title('Bode Plot - Magnitude Response');
xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)');
legend('Location', 'Best');

subplot(2, 1, 2);
grid on;
title('Bode Plot - Phase Response');
xlabel('Frequency (Hz)'); ylabel('Phase (degrees)');
legend('Location', 'Best');

sgtitle('Bode Plots for Different Values of c');

% Polar plot
figure;
for i = 1:length(c_values)
    c = c_values(i);
    G = tf([b, 1], [c, 1, 0, 0]);
    
    % Calculate the Bode plot
    [m, p] = bode(G);
   
    % Make into 1D array
    m = m(1, :);
    p = p(1, :);
    
    % Generate the polar plot
    polarplot((p* pi)/180, m, 'DisplayName', ['c = ' num2str(c)], 'LineWidth', 2);
    hold on;
end
legend('Location', 'Best');
pax = gca;
rlim([0 50]); % set axis limits

%% D //----------//----------//----------//----------//----------//

a = 1;  b = 2;  wn = b;  zeta_values = [0.2, 0.5, 0.8];

% Magnitude and phase plots
figure;
for i = 1:length(zeta_values)
    zeta = zeta_values(i);
    G = tf((a * wn^2), [1, 2 * zeta * wn, wn^2])
    
    % Calculate the Bode plot
    [m, p, w] = bode(G);
    
    % Convert angular frequency to regular frequency (Hz)
    f = w/(2*pi);
    
    % Plot magnitude response
    subplot(2, 1, 1);
    semilogx(f, 20*log10(squeeze(m)), 'DisplayName', ['ζ = ' num2str(zeta)], 'LineWidth', 2);
    hold on;
    
    % Plot phase response
    subplot(2, 1, 2);
    semilogx(f, squeeze(p), 'DisplayName', ['ζ = ' num2str(zeta)], 'LineWidth', 2);
    hold on;
end

% Set labels and legends and show plot
subplot(2, 1, 1);
grid on;
title('Bode Plot - Magnitude Response');
xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)');
legend('Location', 'Best');

subplot(2, 1, 2);
grid on;
title('Bode Plot - Phase Response');
xlabel('Frequency (Hz)'); ylabel('Phase (degrees)');
legend('Location', 'Best');

sgtitle('Bode Plots for Different Values of ζ');

% Polar plot
figure;
for i = 1:length(zeta_values)
    zeta = zeta_values(i);
    G = tf((a * wn^2), [1, 2 * zeta * wn, wn^2]);
    
    % Calculate the Bode plot
    [m, p] = bode(G);
        
    % Make into 1D array
    m = m(1, :);
    p = p(1, :);
    
    % Generate the polar plot
    polarplot((p*pi)/180, m, 'DisplayName', ['ζ = ' num2str(zeta)], 'LineWidth', 2);
    hold on;
end
legend('Location', 'Best');
pax = gca; % set axis limits
rlim([0 3]);
