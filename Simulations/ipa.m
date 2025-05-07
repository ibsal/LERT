function state = ipa(P, T)
    A = 4.57795;
    B = 1221.423;
    C = -87.474;
    vaporPressure = 10^(A - B/(T+C));
    if (100000*vaporPressure)>P
        state = "VAPOR";
    else
        state = "LIQUID";
    end
end
