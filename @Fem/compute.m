function Fem = compute(Fem,varargin)

    p = inputParser;
    addOptional(p,'full',true);
    parse(p,varargin{:});

    if p.Results.full
        Fem = assembleGlobalFem(Fem);
    else
        Fem = assemblePartialFem(Fem);
    end
    
    Fem = assembleBoundaryFem(Fem);
    
    Fem.system.fResidual = Fem.system.fElastic + Fem.system.fDamping ... 
        - Fem.system.fInput - Fem.system.fBody - Fem.system.fContact ...
        - Fem.system.fDilation;
end