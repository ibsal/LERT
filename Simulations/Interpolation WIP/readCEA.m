%=====================================================================
%  readCEA_chamber_v9.m  –  strips stray quote from transport numbers
%---------------------------------------------------------------------
%  Fix: Some CEA tables prepend a single‑quote to millipoise or other
%  transport entries (e.g. '1.2345-4).  toNum() now removes any leading
%  non‑numeric char before parsing so mu_c no longer contains "'".
%---------------------------------------------------------------------
function cea = readCEA(fname)
    raw = fileread(fname);
    blkIdx = regexp(raw,'(?m)^\s*Pin =','start');
    blkIdx(end+1) = numel(raw)+1;

    cols = {'O_F','Pc','T_c','rho_c','Cp_c','gamma_c','a_c', ...
            'mu_c','k_c','Pr_c','Cstar','CF','Isp'};
    for c = cols, data.(c{1}) = []; end

    for p = 1:numel(blkIdx)-1
        Pc = toNum(regexp(raw(blkIdx(p):blkIdx(p)+120),'Pin =\s*([0-9.\'']+)','tokens','once'));
        sec = raw(blkIdx(p):blkIdx(p+1)-1);
        ofIdx = regexp(sec,'(?m)^\s*O/F=','start'); ofIdx(end+1)=numel(sec)+1;
        for k = 1:numel(ofIdx)-1
            sub = sec(ofIdx(k):ofIdx(k+1)-1);
            rec.O_F = toNum(regexp(sub,'O/F=\s*([0-9.\'']+)','tokens','once'));
            rec.Pc = Pc*6894.76;
            rec.T_c   = grab(sub,'(?m)^\s*T,\s*K\s+([0-9.\-\'']+)');
            rec.rho_c = grab(sub,'(?m)^\s*RHO,\s*KG/CU M\s+([0-9.\-\'']+)');
            rec.gamma_c = grab(sub,'(?m)^\s*GAMMAs\s+([0-9.\-\'']+)');
            rec.a_c   = grab(sub,'(?m)^\s*SON VEL,?M/SEC\s+([0-9.\-\'']+)');
            rec.mu_c = grab(sub,'(?m)^\s*VISC,?MILLIPOISE\s+['' ]*([0-9.\-]+)')*1e-4;
           

            eqBlock = regexp(sub,'WITH EQUILIBRIUM REACTIONS[\s\S]*?WITH FROZEN','match','once');
            rec.Cp_c = grab(eqBlock,'Cp, KJ/\(KG\)\(K\)\s+([0-9.\-\'']+)')*1e3;
            rec.k_c  = grab(eqBlock,'(?m)^\s*CONDUCTIVITY\s+['' ]*([0-9.\-]+)')*0.1;
            rec.Pr_c = grab(eqBlock,'PRANDTL NUMBER\s+([0-9.\-\'']+)');

            perf = regexp(sub,'PERFORMANCE PARAMETERS[\s\S]*?$','match','once');
            rec.Cstar = grab(perf,'(?m)^\s*CSTAR,\s*M/SEC\s+([0-9.\-\'']+)');
            rec.CF    = grab(perf,'(?m)^\s*CF\s+([0-9.\-\'']+)');
            rec.Isp   = grab(perf,'(?m)^\s*Isp[^\n]*?\s+([0-9.\-\'']+)');

            for c = cols, data.(c{1})(end+1,1)=rec.(c{1}); end
        end
    end
    cea = struct(); for c=cols, cea.(c{1})=data.(c{1}); end
end

%---------------- helper wrappers -------------------------------------
function num = grab(str,rgx)
    if isempty(str), num = NaN; return; end
    tok = regexp(str,rgx,'tokens','once');
    if isempty(tok), num = NaN; else, num = toNum(tok); end
end

function val = toNum(tok)
    if isempty(tok), val = NaN; return; end
    if iscell(tok), tok = tok{1}; end
    tok = regexprep(tok,'^[^0-9+\-\.]+','');   % strip leading junk
    tok = strtrim(tok);
    % remove leading non‑numeric chars like ' or "+"
    while ~isempty(tok) && ~ismember(tok(1),'-+.0123456789')
        tok = tok(2:end);
    end
    if isempty(tok), val = NaN; return; end
    v = str2double(tok);
    if ~isnan(v), val = v; return; end
    m = regexp(tok,'^([0-9.]+)([\-+])([0-9]+)$','tokens','once');
    if ~isempty(m)
        val = str2double(m{1}) * 10^(str2double([m{2} m{3}])); return; end
    m = regexp(tok,'^([0-9.]+)\s+([0-9]+)$','tokens','once');
    if ~isempty(m)
        val = str2double(m{1}) * 10^(str2double(m{2})); return; end
    val = NaN;
end
