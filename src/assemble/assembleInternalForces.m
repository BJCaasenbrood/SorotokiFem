function Fem = assembleInternalForces(Fem)

beta    = Fem.options.loadingFactor;
[E,~,V] = materialFieldFem(Fem);    
qa      = Fem.system.Ia;

% global matrix
Fem.system.Mass = sparse(Fem.triplets.i,Fem.triplets.j,Fem.triplets.m);
Fem.system.Mass = Fem.system.Mass(qa,qa);

% global damping
Fem.system.Damping = sparse(Fem.triplets.i,Fem.triplets.j,Fem.triplets.c);  
Fem.system.Damping = Fem.system.Damping(qa,qa);
 % init damping force vector
Fem.system.fDamping = Fem.system.Damping * Fem.solver.sol.dx(qa);

% global linear stiffness
Fem.system.Stiffness = sparse(Fem.triplets.i,Fem.triplets.j,           ...
    E(Fem.triplets.e).*Fem.triplets.k);  
Fem.system.Stiffness = Fem.system.Stiffness(qa,qa);

 % global tangent stiffness
Fem.system.Tangent   = sparse(Fem.triplets.i,Fem.triplets.j,           ...
    E(Fem.triplets.e).*Fem.triplets.t); 
Fem.system.Tangent = Fem.system.Tangent(qa,qa);

% global elastic force 
Fem.system.fElastic  = sparse(Fem.triplets.i,1,E(Fem.triplets.e)       ...
    .*Fem.triplets.fi);  

if Fem.options.isNonlinear    
    Fem.system.fElastic = Fem.system.fElastic(qa);
else
    Fem.system.fElastic = Fem.system.Tangent * Fem.solver.sol.x(qa);
end

 % global body force
Fem.system.fBody = sparse(Fem.triplets.i,1,E(Fem.triplets.e)       ...
    .*Fem.triplets.fb);    
Fem.system.fBody = beta * Fem.system.fBody(qa);
    
end
