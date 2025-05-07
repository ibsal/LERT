function v = ipaKinV(T)
    % ipaKinV  Kinematic viscosity of isopropanol at temperature T (K)
    %          Uses a cached curve‑fit so the fit is built only once.
    
    persistent fitFunc                         % cache across calls
    if isempty(fitFunc)                        % build on first call
        tb      = readtable("Kinematic_Viscosity_of_Propan-2-ol_(pure).csv");
        visc    = tb.KinematicViscosity_mm2S_1_;
        temp    = tb.SystemTemperature_K_;
        fitFunc = fit(temp, visc, 'power2');   % create and store the fit
    end
    
    v = fitFunc(T) * 10^-6;                            % evaluate the cached fit
end
