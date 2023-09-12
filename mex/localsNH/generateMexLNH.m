cfg = coder.config('mex');
cfg.IntegrityChecks = false;
cfg.ExtrinsicCalls = false;
cfg.IntegrityChecks = false;
cfg.SaturateOnIntegerOverflow = false;
cfg.ResponsivenessChecks = false;
cfg.NumberOfCpuThreads = 16;

codegen -config cfg LocalsNHFast ...
    -args {coder.typeof(10,[inf,inf]), coder.typeof(10,[inf,inf]), coder.typeof(10,[1,1]), coder.typeof(10,[1,1]), coder.typeof(10,[1,1]), coder.typeof(10,[inf,inf]), coder.typeof(10,[inf,1,inf]), coder.typeof(10,[inf,inf,inf]), coder.typeof(10,[inf,1]), coder.typeof(10,[inf,inf]), coder.typeof(10,[inf,inf]), coder.typeof(10,[1,1]), coder.typeof(10,[1,1]), coder.typeof(10,[inf,1]), coder.typeof(10,[1,1]), coder.typeof(10,[1,1])}

codegen -config cfg LocalsNHFastElastic -config:mex ...
    -args {coder.typeof(10,[inf,inf]), coder.typeof(10,[inf,inf]), coder.typeof(10,[1,1]), coder.typeof(10,[1,1]), coder.typeof(10,[1,1]), coder.typeof(10,[inf,inf]), coder.typeof(10,[inf,1,inf]), coder.typeof(10,[inf,inf,inf]), coder.typeof(10,[inf,1]), coder.typeof(10,[inf,inf]), coder.typeof(10,[inf,inf]), coder.typeof(10,[1,1]), coder.typeof(10,[1,1]), coder.typeof(10,[inf,1]), coder.typeof(10,[1,1]), coder.typeof(10,[1,1])}