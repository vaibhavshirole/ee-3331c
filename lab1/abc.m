%% 2.1

% G = 1/(5*s+2) * e^(-2*s)

A = 1;
tau = 5;
td = 2;
Tend = 50;

Gp = tf(A,[tau 2],'InputDelay', td);

[y,t] = step(Gp, Tend);
figure(1), plot(t,y),grid,title('Step response of plant');
xlabel('Time (s)');
ylabel('Amplitude');

% Add labels for A, td, and tau
text(1, 1/2, ['A = ', num2str(1/2)], 'FontSize', 12, 'VerticalAlignment', 'top');
text(td, 0, ['td = ', num2str(td)], 'FontSize', 12, 'VerticalAlignment', 'top');
text(tau, 1/2*(63/100), ['tau = ', num2str(tau)], 'FontSize', 12, 'VerticalAlignment', 'bottom');

%% 2.2 - P control
x = 1/((A/tau)*td);
h = 1.25;
l = 0.75;

Kp = x;
Gol = series(Kp, Gp);
Gcl = feedback(Gol, 1);
[y1, t1] = step(Gcl, Tend);

Kp = x*h;
Gol = series(Kp, Gp);
Gcl = feedback(Gol, 1);
[y2, t2] = step(Gcl, Tend);

Kp = x*l;
Gol = series(Kp, Gp);
Gcl = feedback(Gol, 1);
[y3, t3] = step(Gcl, Tend);

figure(2),plot(t1,y1,t2,y2,t3,y3), grid;
title('Closed loop step response with P-control');
xlabel('Time (s)');
ylabel('Amplitude');
legend('x', 'xh', 'xl');


[mag,php,w] = bode(Gp);
dbp=20*log10(mag);

Kp = x;
[mag,phol1] = bode(series(Kp,Gp),w);
dbol1=20*log10(mag);

Kp = x*h;
[mag,phol2] = bode(series(Kp,Gp),w);
dbol2 = 20*log10(mag);

figure(3);
subplot(211), semilogx(w, dbp(:), w, dbol1(:), w, dbol2(:)), grid;
title('Open loop Bode (mag) plot with P−Control');
legend('Plant','gain=Kp','gain>Kp');
xlabel('Frequency (rad/s)');
ylabel('Magnitude (dB)');
subplot(212), semilogx(w, php(:), w, phol1(:), w, phol2(:)), grid;
title('Open loop Bode (phase) plot with P−Control');
legend('Plant','gain=Kp','gain>Kp');
xlabel('Frequency (rad/s)');
ylabel('Phase (degrees)');


Kp = x;
[mag,phcl1] = bode(feedback(series(Kp,Gp),1),w);
dbcl1=20*log10(mag);

Kp = x*h;
[mag,phcl2] = bode(feedback(series(Kp,Gp),1),w);
dbcl2=20*log10(mag);

figure(4);
subplot(211), semilogx(w, dbcl1(:), w, dbcl2(:)), grid;
title('Closed loop Bode (mag) plot with P−Control');
legend('gain=Kp','gain>Kp');
xlabel('Frequency (rad/s)');
ylabel('Magnitude (dB)');
subplot(212), semilogx(w, phcl1(:), w, phcl2(:)), grid;
title('Closed loop Bode (phase) plot with P−Control');
legend('gain=Kp','gain>Kp');
xlabel('Frequency (rad/s)');
ylabel('Phase (degrees)');

%% 2.3 - PI control
alpha = 0.9/((A/tau)*td);
beta = td/0.3;

Kp = alpha;
Ti = beta;
C = Kp*tf([Ti 1],[Ti 0]);
Gol = series(C,Gp);
Gcl = feedback(Gol,1);
[y4,t4] = step(Gcl, Tend);

Kp = alpha*h;
Ti = beta;
C = Kp*tf([Ti 1],[Ti 0]);
Gol = series(C,Gp);
Gcl = feedback(Gol,1);
[y5,t5] = step(Gcl, Tend);

figure(5),plot(t4,y4,t5,y5),grid;
title('Closed loop step response with PI-control');
legend('gain=Kp','gain>Kp');
xlabel('Time (s)');
ylabel('Amplitude');


[mag,php,w] = bode(Gp);
dbp=20*log10(mag);

Kp = alpha;
Ti = beta;
C = Kp*tf([Ti 1],[Ti 0]);
[mag, phol1] = bode(series(C, Gp), w);
dbol1 = 20 * log10(mag);

Kp = alpha * h;
Ti = beta;
C = Kp*tf([Ti 1],[Ti 0]);
[mag, phol2] = bode(series(C, Gp), w);
dbol2 = 20 * log10(mag);

figure(6);
subplot(211), semilogx(w, dbp(:), w, dbol1(:)), grid;
title('Open-loop Bode (mag) plot with PI-Control');
legend('Plant', 'gain=Kp');
xlabel('Frequency (rad/s)');
ylabel('Magnitude (dB)');
subplot(212), semilogx(w, php(:), w, phol1(:)), grid;
title('Open-loop Bode (phase) plot with PI-Control');
legend('Plant', 'gain=Kp');
xlabel('Frequency (rad/s)');
ylabel('Phase (degrees)');