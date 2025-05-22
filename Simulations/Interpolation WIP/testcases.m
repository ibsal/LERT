%% quickSanity.m  — verifies parser + grid interpolator
%  Assumes:
%     • rocketCEAgrid.csv       already exists in pwd
%     • ceaGridInterp.m         on path
%

T = readtable('rocketCEAgrid.csv');     % full grid

%% ---------- 1.  exact-value checks for first Pin/O-F block -------------
% Pin 14 psia  → Pc_bar = 0.96526 ;  O/F = 1.0  (see txt file)
tol = 1e-1;

row = T( abs(T.Pc-96526.64) < 1e-1  &  abs(T.O_F-1.0) < tol , : );

assert(~isempty(row),'Row O/F=1, Pc=0.96526 bar not found');

tests = { ...
    'T_c',    1351.17 ;                                    % K
    'rho_c',  0.13393 ;                                    % kg/m^3   1.3393-1
    'Cp_c',   2435.5 ;                                     % J/(kg·K) 2.4355 kJ
    'gamma_c',1.2997 ;
    'a_c',    967.9 ;                                      % m/s
    'mu_c',   0.70875e-4 ;                                 % Pa·s
    'k_c',    0.53615 ;                                    % W/(m·K)  5.3615 mW/cm·K
    'Pr_c',   0.3220 ;
    'Cstar',  1275.2 ;                                     % m/s
    'CF',     0.7038 ;
    'Isp',    897.4                                        % m/s  (CEA SI)
};

for k = 1:size(tests,1)
    fld = tests{k,1};
    ref = tests{k,2};
    val = row.(fld);
    assert( abs(val-ref) < 1e-3*abs(ref)+1e-9 , ...
        'Mismatch in %s  (got %.5g, expected %.5g)',fld,val,ref);
end
disp('✓  Chamber block @14 psia, O/F=1  — all values OK');

%% ---------- 2.  interpolation spot checks ------------------------------
% Example:  Pc = 150 bar   O/F = 3.5
Pc  = 150000;   OF = 3.5;

Tc_interp = ceaGridInterp('rocketCEAgrid.csv','T_c',  Pc, OF);
mu_interp = ceaGridInterp('rocketCEAgrid.csv','mu_c', Pc, OF);

fprintf('Interpolated T_c(%.0f bar, O/F=%.1f)  = %.2f K\\n',Pc,OF,Tc_interp);
fprintf('Interpolated mu_c(%.0f bar, O/F=%.1f) = %.3e Pa·s\\n',Pc,OF,mu_interp);

% quick reasonableness:
assert(Tc_interp>1500 & Tc_interp<4000,'T_c interpolation out of range');
assert(mu_interp>5e-5  & mu_interp<2e-4,'mu_c interpolation out of range');

disp('✓  Interpolation spot check passed');
