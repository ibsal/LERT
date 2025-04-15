v = @(M) sqrt((gamma+1)/(gamma-1)) * atan(sqrt(((gamma-1)/(gamma+1)) * (M^2-1))) - atan(sqrt(M^2-1));

M1 = 2.6;

alpha = 15;

gamma = 1.4;

zerome = @(M) rad2deg(v(M)-v(M1)) - alpha;

M2 = fzero(zerome, M1);

P2divP1 = (((1 + 0.5*(gamma-1)*M2^2)/(1 + 0.5*(gamma-1)*M1^2))^(-gamma/(gamma-1)));

%annonymous function for beta theta m equation
bsolv = @(beta, theta, M1, y) ((2 * cot(beta) * (M1.^2 * sin(beta).^2 - 1))./ ...
 (M1.^2 * (y + cos(2*beta)) + 2)) - tan(theta);
%parameters
theta = alpha;
%anonymous function for fzero using parameters
fun = @(x) bsolv(x, theta * pi/180, M1, gamma);
%strong shock with initial guess at 90 degrees
bstrong = 180 * fzero(fun, 90 * 0.01745329251)/pi;
%weak shock between 0 and strong shock angle
bweak = 180 * fzero(fun, [0.01, bstrong*0.01745329251])/pi;
%Assume weak shocks:
b = bweak;
Mn1 = M1 * sin(bweak*0.01745329251);
Mn3 = sqrt((1 + Mn1^2 * 0.5 * (gamma-1))/(gamma * Mn1^2 - 0.5*(gamma-1)));
rhoratio = ((gamma+1)*Mn1^2 )/ (2+ (gamma-1) * Mn1^2);
pratio = 1 + ((2*gamma)/(gamma + 1)) * (Mn1^2 - 1);
M3 = Mn3/sin((b-theta)*0.01745329251);
pstagratio = pratio * ((1 + 0.5*(gamma-1) * Mn3^2)/(1 + 0.5*(gamma-1) * Mn3^2))^(gamma/(gamma-1));
tratio = (1 + 0.5*(gamma-1) * M1^2)/(1 + 0.5*(gamma-1) * M3^2);

Cl = (2/(gamma * M1^2)) * (pratio - P2divP1) * cos(deg2rad(alpha))
Cd = (2/(gamma * M1^2)) * (pratio - P2divP1) * sin(deg2rad(alpha))
