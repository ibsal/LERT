PT = 300;
P = 159;
gamma = 1.4;
mexit = @(M) (1 + 0.5*(gamma-1)*M^2)^((-gamma)/(gamma-1)) -  P/PT;
pcrit = @(crit) (1 + 0.5*(gamma-1)*1^2)^((-gamma)/(gamma-1)) -  crit/PT;
try
    fzero(mexit, [1e-6, 0.9999999999999])
end
try
    fzero(mexit, [1.00000000000001, 100])
end

fzero(pcrit, P)

A1/A* = 