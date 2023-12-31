% A0271121X
% a = 1
%
% G(s) = 1/(s(s+a)), K(s) = K
% velocity error constant, Kv = 20
% phase margin of at least 50 deg

%% PART A - Lead compensator (see work in assignment)

D1 = tf([0.34, 1], [0.066, 1]); % +5 deg added to get >50 deg target

%% PART B - Design a lag compensator to meet the above specifications

D2 = tf([17.33, 1], [519.93, 1]); % +10 deg added to get >50 deg target

%% PART C - Simulate your results in Matlab and show the Bode plots

% Uncompensated system G(s)*K(s)
GK = tf(20, [1, 1, 0]); % numerator, denominator
figure;
h1 = bodeplot(GK);
grid on;

% Gain and phase margin analysis for G(s)*K(s)
[Gm, Pm, Wcg, Wcp] = margin(GK);
title(['{NONE}  Gm: ', num2str(Gm), ' dB (at ',num2str(Wcg),' rad/s), Pm: ', ...
    num2str(Pm), ' deg (at ',num2str(Wcp),' rad/s)']);


% Lead compensated system G(s)*D1(s)*K(s)
GD1K = GK*D1;
figure;
h2 = bodeplot(GD1K);
grid on;

% Gain and phase margin analysis for G(s)*D1(s)*K(s)
[Gm, Pm, Wcg, Wcp] = margin(GD1K);
title(['{LEAD}  Gm: ', num2str(Gm), ' dB (at ',num2str(Wcg),' rad/s), Pm: ', ...
    num2str(Pm), ' deg (at ',num2str(Wcp),' rad/s)']);


% Lag compensated system G(s)*D2(s)*K(s)
GD2K = GK*D2;
figure;
h3 = bodeplot(GD2K);
grid on;

% Gain and phase margin analysis for G(s)*D2(s)*K(s)
[Gm, Pm, Wcg, Wcp] = margin(GD2K);
title(['{LAG}  Gm: ', num2str(Gm), ' dB (at ',num2str(Wcg),' rad/s), Pm: ', ...
    num2str(Pm), ' deg (at ',num2str(Wcp),' rad/s)']);
