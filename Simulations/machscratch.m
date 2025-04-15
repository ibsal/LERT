P2 = 0.1306;
P1 = 1;
M1 = 1.58;
gamma = 1.4;


fun = @(M2) ((1 + 0.5*(gamma-1)*M2^2)/(1 + 0.5*(gamma-1)*M1^2))^(-gamma/(gamma-1)) - P2/P1;

mtwo = fzero(fun, 0.5);

fprintf("%.6f\n", mtwo)

v = @(M) sqrt((gamma+1)/(gamma-1)) * atan(sqrt(((gamma-1)/(gamma+1)) * (M^2-1))) - atan(sqrt(M^2-1));

theta = rad2deg(v(mtwo)-v(M1));

fprintf("%.6f\n", theta)

%annonymous function for beta theta m equation
bsolv = @(beta, theta, M1, y) ((2 * cot(beta) * (M1.^2 * sin(beta).^2 - 1))./ ...
 (M1.^2 * (y + cos(2*beta)) + 2)) - tan(theta);
%parameters
theta = 30.6;
M1 = 3;
gamma = 1.4;
%anonymous function for fzero using parameters
fun = @(x) bsolv(x, theta * pi/180, M1, gamma);
%strong shock with initial guess at 90 degrees
bstrong = 180 * fzero(fun, 90 * 0.01745329251)/pi;
%weak shock between 0 and strong shock angle
bweak = 180 * fzero(fun, [0.01, bstrong*0.01745329251])/pi;
%Assume weak shocks:
b = bweak;
Mn1 = M1 * sin(bweak*0.01745329251);
Mn2 = sqrt((1 + Mn1^2 * 0.5 * (gamma-1))/(gamma * Mn1^2 - 0.5*(gamma-1)));
rhoratio = ((gamma+1)*Mn1^2 )/ (2+ (gamma-1) * Mn1^2);
pratio = 1 + ((2*gamma)/(gamma + 1)) * (Mn1^2 - 1);
M2 = Mn2/sin((b-theta)*0.01745329251);
pstagratio = pratio * ((1 + 0.5*(gamma-1) * Mn2^2)/(1 + 0.5*(gamma-1) * Mn1^2))^(gamma/(gamma-1));
tratio = (1 + 0.5*(gamma-1) * M1^2)/(1 + 0.5*(gamma-1) * M2^2);
%make results nice
fprintf("\nParameters:\nM1 = %.6f, Turn Angle = %.6f, Gamma = %.6f\n\n", M1, theta, gamma);
fprintf("Results:\nBeta = %.6f, Mn1 = %.6f, Mn2 = %.6f, M2 = %.6f\n" + ...
 "Rho2/Rho1 = %.6f, P2/P1 = %.6f, Po2/Po1 = %.6f, T2/T1 = %.6f\n\n", ...
 b, Mn1, Mn2, M2, rhoratio, pratio, pstagratio, tratio);
%NORMAL SHOCKS
M2 = sqrt((1 + M1^2 * 0.5 * (gamma-1))/(gamma * M1^2 - 0.5*(gamma-1)));
pratio = 1 + ((2*gamma)/(gamma + 1)) * (M1^2 - 1);
pstagratio = pratio * ((1 + 0.5*(gamma-1) * M2^2)/(1 + 0.5*(gamma-1) * M1^2))^(gamma/(gamma-1));
fprintf("\n\nNormal Shock:\nM2 = %.6f, Po2/Po1 = %.6f, P2/P1 = %.6f\n\n", M2, pstagratio, pratio)

zerome = @(x) rad2deg(v(x) - v(1.363372)) - 30.6;

m3 = fzero(zerome, 2);

fprintf("%.6f\n", m3)

fprintf("%.6f\n", ((1 + 0.5*(gamma-1)*m3^2)/(1 + 0.5*(gamma-1)*1.363372^2))^(-1))