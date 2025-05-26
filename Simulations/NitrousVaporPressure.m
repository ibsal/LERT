function P = NitrousVaporPressure(T)
    A = 9.3867;
    B =  748.71;
    C = -13.306;
    P = (10^(A - (B/(T+C))));

end
