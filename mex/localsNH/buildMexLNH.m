function flag = buildMexLNH
    flag = 0;
    try
    cfg = coder.config('mex');
    cfg.IntegrityChecks = false;
    cfg.ExtrinsicCalls = false;
    cfg.IntegrityChecks = false;
    cfg.SaturateOnIntegerOverflow = false;
    cfg.ResponsivenessChecks = false;
    cfg.NumberOfCpuThreads = 16;

    warning off
    codegen -config cfg LocalsNHFast -silent ...
        -args {coder.typeof(10,[inf,inf]), coder.typeof(10,[inf,inf]), coder.typeof(10,[1,1]), coder.typeof(10,[1,1]), coder.typeof(10,[1,1]), coder.typeof(10,[inf,inf]), coder.typeof(10,[inf,1,inf]), coder.typeof(10,[inf,inf,inf]), coder.typeof(10,[inf,1]), coder.typeof(10,[inf,inf]), coder.typeof(10,[inf,inf]), coder.typeof(10,[1,1]), coder.typeof(10,[1,1]), coder.typeof(10,[inf,1]), coder.typeof(10,[1,1]), coder.typeof(10,[1,1])}
        
    warning on    
    flag = 1;
    end
end