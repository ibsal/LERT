%% CEA block‑parser & interpolator (v4 – Cstar, Mach_e, a_e) -------------------
%  readCEA   – now also extracts exit Mach number and exit sonic velocity
%  ceaInterp – unchanged (scalar return fix retained)
%  ------------------------------------------------------------------
%  New struct fields returned by readCEA:
%      • Mach_e   – exit Mach number (dimensionless)
%      • a_e      – sonic velocity at exit [m/s]
% -------------------------------------------------------------------

function cea = readCEA(fname)
    arguments, fname (1,:) char, end
    txt = fileread(fname);

    blkIdx = regexp(txt,'(?m)^\s*O/F=','start');
    if isempty(blkIdx)
        error("readCEA:NoBlocks","No 'O/F=' markers in %s", fname);
    end
    n = numel(blkIdx);

    vars = {'O_F','Pc','Tc','Te','gamma_e','AeAt','CF','Isp','Cstar','Mach_e','a_e'};
    for v = vars, data.(v{1}) = nan(n,1); end

    blkIdx(end+1) = numel(txt)+1;   % sentinel

    for k = 1:n
        d = readOneBlock( txt( blkIdx(k):blkIdx(k+1)-1 ) );
        for v = vars, data.(v{1})(k) = d.(v{1}); end
    end

    cea = struct('O_F',data.O_F,      'Pc_bar',data.Pc,   'T_c',data.Tc,  ...
                 'T_e',data.Te,       'gamma_e',data.gamma_e, 'Ae_At',data.AeAt, ...
                 'CF',data.CF,        'Isp',data.Isp,     'Cstar',data.Cstar, ...
                 'Mach_e',data.Mach_e,'a_e',data.a_e);
end

% ---------------- block parser -------------------------------------------
function blk = readOneBlock(str)
    blk.O_F   = grab(str,'O/F=\s*([0-9.]+)');
    blk.Pc    = grab(str,'(?m)^\s*P,\s*BAR\s+([0-9.eE+-]+)');

    t = regexp(str,'(?m)^\s*T,\s*K\s+([0-9.eE+-]+)\s+[0-9.eE+-]+\s+([0-9.eE+-]+)','tokens','once');
    blk.Tc = str2double(t{1});  blk.Te = str2double(t{2});

    blk.gamma_e = grab(str,'(?m)^\s*GAMMAs\s+[0-9.eE+-]+\s+[0-9.eE+-]+\s+([0-9.eE+-]+)');
    blk.AeAt    = grab(str,'(?m)^\s*Ae/At\s+[0-9.eE+-]+\s+([0-9.eE+-]+)');
    blk.CF      = grab(str,'(?m)^\s*CF\s+[0-9.eE+-]+\s+([0-9.eE+-]+)');
    blk.Isp     = grab(str,'(?m)^\s*Isp[^\n]*?\s+[0-9.eE+-]+\s+([0-9.eE+-]+)');
    blk.Cstar   = grab(str,'(?m)^\s*CSTAR,\s*M/SEC\s+([0-9.eE+-]+)');

    blk.a_e     = grab(str,'(?m)^\s*SON VEL,?M/SEC\s+[0-9.eE+-]+\s+[0-9.eE+-]+\s+([0-9.eE+-]+)');
    blk.Mach_e  = grab(str,'(?m)^\s*MACH NUMBER\s+[0-9.eE+-]+\s+[0-9.eE+-]+\s+([0-9.eE+-]+)');
end

% ----- helper ------------------------------------------------------------
function val = grab(str, rgx)
    tok = regexp(str, rgx, 'tokens', 'once');
    if isempty(tok)
        val = NaN;   % keep NaN if not found
    else
        val = str2double(tok{1});
    end
end