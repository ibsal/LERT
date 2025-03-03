x = 0.1016
gamma = 0.1400
RadiusAtX = 0.0508
Rthroat = 0.0127


fun = @(M) (1/M^2) * ((2/(gamma + 1)) * (1 + 0.5*(gamma-1)*M^2))^((gamma + 1)/(gamma-1)) - (RadiusAtX/Rthroat)^2;

Mach = fzero(fun, 2);