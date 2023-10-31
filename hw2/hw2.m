%% A //----------//----------//----------//----------//----------//

a = 1;
b = 2;
b_values = [0.05*b, 1*b, 50*b];

figure;
hold on;
for i = 1:length(b_values)
    b = b_values(i);
    G = tf(a, [1, b])
    bode(G);
end
hold off;

grid on;
legend('0.05b', 'b', '50b');
title('Bode Plots for Different Values of b');

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
hold on;
for i = 1:length(c_values)
    c = c_values(i);
    G = tf([b, 1], [c, 1, 0, 0])
    bode(G);
end
hold off;

grid on;
legend('0.1b', '10b');
title('Bode Plots for Different Values of c');

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
hold on;
for i = 1:length(zeta_values)
    zeta = zeta_values(i);
    G = tf((a * wn^2), [1, 2 * zeta * wn, wn^2])
    bode(G);
end
hold off;

grid on;
legend('0.2', '0.5', '0.8');
title('Bode Plots for Different Values of ζ');

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
