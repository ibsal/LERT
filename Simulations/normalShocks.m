%% Normal Shock Solver
M2 = @(M) sqrt((1 + ((gamma-1)/2) * M^2)/(gamma * M^2 - (gamma -1)/2));
Pratio = @(M) 1 + ((2*gamma)/(gamma+1)) * (M^2 - 1);