function f = darcy(e, D, Re)
    func = @(f) 1/sqrt(f) + 2*log10((e/D)/3.7 + 2.51/(Re * sqrt(f)));
    f = fzero(func, [1e-7, 100]);
end