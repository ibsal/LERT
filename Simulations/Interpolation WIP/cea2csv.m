%=====================================================================
%  CEA CSV pipeline  (v2) – *regex‑free* full‑field exporter
%---------------------------------------------------------------------
%  Changes in v2 (2025‑05‑21)
%    •  No manual regex needed: uses the robust readCEA() parser to
%       extract *every* performance / transport field already defined
%       there (O_F, Pc_bar, T_c, T_e, gamma_c, gamma_e, Ae_At, CF, Isp,
%       Cstar, Mach_e, a_e, Cp, Pr, mu_e … add more in readCEA once and
%       they flow into the CSV automatically).
%    •  cea2csv now one‑liners: T = struct2table(readCEA(..));
%    •  ceaGridInterp unchanged – it loads whatever columns exist.
%---------------------------------------------------------------------
%  Usage:
%       cea2csv('rocketceacomplete.txt','rocketCEAgrid.csv');
%       Tc = ceaGridInterp('rocketCEAgrid.csv','T_c', Pc_bar, OF);
%=====================================================================

function cea2csv(infile,outfile)
% Parse CEA text → CSV using readCEA (no regex here)
if nargin<2, outfile='rocketCEAgrid.csv'; end
S = readCEA(infile);          % <-- your v6 parser with full field list
T = struct2table(S);
writetable(T,outfile);
fprintf('Created %s  (%d rows, %d columns)\n',outfile,height(T),width(T));
end

%---------------------------------------------------------------------
