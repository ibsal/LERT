%=========================================================================
% readCEA.m  (v6)  –  parse NASA CEA short‑format tables
%-------------------------------------------------------------------------
%  NEW in v6
%    •  Adds *combustion‑zone γ*  (gamma_c) to the output struct.
%    •  Keeps the previously parsed *exit γ* (gamma_e).
%    •  Cp and Pr exit transport terms retained from v5.
%-------------------------------------------------------------------------
%  Usage:
%       cea = readCEA('rocketN2OxIPA.txt');
%       Tc  = ceaInterp(cea,'T_c','O_F',3.0);
%       gc  = ceaInterp(cea,'gamma_c','O_F',3.0);   % <- chamber γ
%-------------------------------------------------------------------------
function cea = readCEA(fname)
    arguments (Input) fname (1,:) char 
    end

    txt   = fileread(fname);
    blkIdx = regexp(txt,'(?m)^\s*O/F=','start');
    if isempty(blkIdx)
        error("readCEA:NoBlocks","No 'O/F=' markers in %s",fname);
    end
    n = numel(blkIdx);
    blkIdx(end+1) = numel(txt)+1;   % sentinel for last block slice

    % ---- pre‑allocate data arrays ------------------------------------
    vars = {'O_F','Pc','Tc','Te','gamma_c','gamma_e','AeAt','CF','Isp', ...
            'Cstar','Mach_e','a_e','Cp','Pr','mu_e'};
    for v = vars, data.(v{1}) = nan(n,1); end

    % ---- loop over each performance block ----------------------------
    for k = 1:n
        sub = txt( blkIdx(k) : blkIdx(k+1)-1 );
        blk = parseBlock(sub);
        for v = vars, data.(v{1})(k) = blk.(v{1}); end
    end

    % ---- assemble output struct -------------------------------------
    cea = struct( ...
        'O_F',      data.O_F , ...
        'Pc_bar',   data.Pc  , ...
        'T_c',      data.Tc  , ...
        'T_e',      data.Te  , ...
        'gamma_c',  data.gamma_c , ...
        'gamma_e',  data.gamma_e , ...
        'Ae_At',    data.AeAt , ...
        'CF',       data.CF , ...
        'Isp',      data.Isp , ...
        'Cstar',    data.Cstar , ...
        'Mach_e',   data.Mach_e , ...
        'a_e',      data.a_e , ...
        'Cp',       data.Cp , ...
        'Pr',       data.Pr , ...
        'mu_e',     data.mu_e );
end

%======================================================================
function blk = parseBlock(str)
% Extract chamber, throat, and exit values from one O/F block

    blk.O_F = grab1(str,'O/F=\s*([0-9.]+)');
    blk.Pc  = grab1(str,'(?m)^\s*P,\s*BAR\s+([0-9.Ee+-]+)');

    t = regexp(str,'(?m)^\s*T,\s*K\s+([0-9.Ee+-]+)\s+[0-9.Ee+-]+\s+([0-9.Ee+-]+)','tokens','once');
    blk.Tc = str2double(t{1});   blk.Te = str2double(t{2});

    % --- GAMMAs line:   chamber  throat  exit -----------------------
    g = regexp(str,'(?m)^\s*GAMMAs\s+([0-9.Ee+-]+)\s+[0-9.Ee+-]+\s+([0-9.Ee+-]+)','tokens','once');
    blk.gamma_c = str2double(g{1});       % chamber γ
    blk.gamma_e = str2double(g{2});       % exit γ  (legacy field)

    blk.AeAt   = grab1(str,'(?m)^\s*Ae/At\s+[0-9.Ee+-]+\s+([0-9.Ee+-]+)');
    blk.CF     = grab1(str,'(?m)^\s*CF\s+[0-9.Ee+-]+\s+([0-9.Ee+-]+)');
    blk.Isp    = grab1(str,'(?m)^\s*Isp[^\n]*?\s+[0-9.Ee+-]+\s+([0-9.Ee+-]+)');
    blk.Cstar  = grab1(str,'(?m)^\s*CSTAR,\s*M/SEC\s+([0-9.Ee+-]+)');
    blk.a_e    = grab1(str,'(?m)^\s*SON VEL,?M/SEC\s+[0-9.Ee+-]+\s+[0-9.Ee+-]+\s+([0-9.Ee+-]+)');
    blk.Mach_e = grab1(str,'(?m)^\s*MACH NUMBER\s+[0-9.Ee+-]+\s+[0-9.Ee+-]+\s+([0-9.Ee+-]+)');

    % ----- transport (exit, equilibrium) ----------------------------
    blk.Cp = grabN(str,'(?m)^\s*Cp,\s*KJ/\(KG\)\(K\)\s+','eq',1) * 1e3;  % kJ→J
    blk.Pr = grabN(str,'(?m)^\s*PRANDTL NUMBER\s+','eq',1);
    mu_mP  = grabN(str,'(?m)^\s*VISC,\s*MILLIPOISE\s+','',1);   % mP → Pa·s
    blk.mu_e = mu_mP * 1e-4;
end

%===================== helper regex wrappers =========================
function val = grab1(str,rgx)
    tok = regexp(str,rgx,'tokens','once');
    val = str2double(tok{1});
end

function num = grabN(str,tag,selector,col)
% Grab the COL‑th numeric after TAG. If selector='eq', start search
% after the line containing that selector keyword (e.g. 'WITH EQUILIBRIUM').
    if ~isempty(selector)
        base = regexp(str,['WITH\s+' upper(selector)],'start','once');
        block = str(base:end);
    else
        block = str;
    end
    pat = [tag '([0-9.Ee+-]+)'];
    for k = 2:col, pat = [pat '\s+([0-9.Ee+-]+)']; end %#ok<AGROW>
    tok = regexp(block,pat,'tokens','once');
    num = isempty(tok) * NaN + ~isempty(tok) * str2double(tok{end});
end
