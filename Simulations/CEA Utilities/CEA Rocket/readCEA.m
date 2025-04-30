%% CEA tools + Mach utilities (v5 – adds Cp & Pr transport) ---------------------
%  readCEA now extracts specific heat (Cp) and Prandtl number (Pr) from
%  each performance block (exit column). New struct fields:
%      • Cp   – exit cp  [J/(kg·K)]
%      • Pr   – exit Prandtl number
%  -----------------------------------------------------------------------

function cea = readCEA(fname)
    arguments (Input) fname (1,:) char
    end
    txt = fileread(fname);
    blkIdx = regexp(txt,'(?m)^\s*O/F=','start');
    if isempty(blkIdx), error("readCEA:NoBlocks","No 'O/F=' markers in %s",fname); end
    n = numel(blkIdx);

    vars = {'O_F','Pc','Tc','Te','gamma_e','AeAt','CF','Isp','Cstar', ...
        'Mach_e','a_e','Cp','Pr','mu_e'};   %  <-- Cp  Pr  (not Cp_e/Pr_e)
    for v = vars, data.(v{1}) = nan(n,1); end
    blkIdx(end+1) = numel(txt)+1;

    for k = 1:n
        blk = parseBlock(txt(blkIdx(k):blkIdx(k+1)-1));
        for v = vars, data.(v{1})(k) = blk.(v{1}); end
    end

    cea = struct('O_F',data.O_F,'Pc_bar',data.Pc,'T_c',data.Tc,'T_e',data.Te,...
             'gamma_e',data.gamma_e,'Ae_At',data.AeAt,'CF',data.CF,...
             'Isp',data.Isp,'Cstar',data.Cstar,'Mach_e',data.Mach_e,...
             'a_e',data.a_e,'Cp',data.Cp,'Pr',data.Pr,'mu_e',data.mu_e);
end

% -------------------------- block parser ---------------------------------
function blk = parseBlock(str)
    blk.O_F   = grab1(str,'O/F=\s*([0-9.]+)');
    blk.Pc    = grab1(str,'(?m)^\s*P,\s*BAR\s+([0-9.Ee+-]+)');

    t = regexp(str,'(?m)^\s*T,\s*K\s+([0-9.Ee+-]+)\s+[0-9.Ee+-]+\s+([0-9.Ee+-]+)','tokens','once');
    blk.Tc = str2double(t{1}); blk.Te = str2double(t{2});

    blk.gamma_e = grab1(str,'(?m)^\s*GAMMAs\s+[0-9.Ee+-]+\s+[0-9.Ee+-]+\s+([0-9.Ee+-]+)');
    blk.AeAt    = grab1(str,'(?m)^\s*Ae/At\s+[0-9.Ee+-]+\s+([0-9.Ee+-]+)');
    blk.CF      = grab1(str,'(?m)^\s*CF\s+[0-9.Ee+-]+\s+([0-9.Ee+-]+)');
    blk.Isp     = grab1(str,'(?m)^\s*Isp[^\n]*?\s+[0-9.Ee+-]+\s+([0-9.Ee+-]+)');
    blk.Cstar   = grab1(str,'(?m)^\s*CSTAR,\s*M/SEC\s+([0-9.Ee+-]+)');
    blk.a_e     = grab1(str,'(?m)^\s*SON VEL,?M/SEC\s+[0-9.Ee+-]+\s+[0-9.Ee+-]+\s+([0-9.Ee+-]+)');
    blk.Mach_e  = grab1(str,'(?m)^\s*MACH NUMBER\s+[0-9.Ee+-]+\s+[0-9.Ee+-]+\s+([0-9.Ee+-]+)');

    % ------------ transport properties (exit column, equilibrium) --------
    blk.Cp  = grabN(str,'(?m)^\s*Cp,\s*KJ/\(KG\)\(K\)\s+','eq',1);   % kJ/kg·K
    blk.Pr  = grabN(str,'(?m)^\s*PRANDTL NUMBER\s+','eq',1);
    mu_mP    = grabN(str,'(?m)^\s*VISC,\s*MILLIPOISE\s+','',1);       % millipoise → Pa·s
    blk.mu_e = mu_mP * 1e-4;                                 % 1 mP = 1e-3 P = 0.1 Pa·s
    blk.Cp  = blk.Cp * 1e3;                                     % kJ → J
end

% ---------- helpers ------------------------------------------------------
function val = grab1(str,rgx), tok=regexp(str,rgx,'tokens','once'); val=str2double(tok{1}); end

function num = grabN(str,tag,selector,col)
% Find a line that starts with TAG. If selector='eq' choose the first
% instance *after* 'WITH EQUILIBRIUM', otherwise just the first.
    if ~isempty(selector)
        base = regexp(str,['WITH\s+' upper(selector)],'start','once');
        block = str(base:end);
    else
        block = str;
    end
    pat = [tag '([0-9.Ee+-]+)'];           % first column
    for k=2:col                            % advance to Nth column
        pat = [pat '\s+([0-9.Ee+-]+)'];
    end
    tok = regexp(block,pat,'tokens','once');
    if isempty(tok), num=NaN; else, num=str2double(tok{end}); end
end

