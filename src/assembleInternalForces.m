function Fem = assembleInternalForces(Fem)

beta    = Fem.options.loadingFactor;
[E,~,V] = materialFieldFem(Fem);    
qa      = Fem.system.Ia;

%M  = sparse(Fem.i,Fem.j,E(Fem.e).*Fem.m);  % init mass matrix
Fem.system.Mass = sparse(Fem.triplets.i,Fem.triplets.j,Fem.triplets.m);
Fem.system.Mass = Fem.system.Mass(qa,qa);

% init mass matrix
Fem.system.Damping = sparse(Fem.triplets.i,Fem.triplets.j,Fem.triplets.c);  
Fem.system.Damping = Fem.system.Damping(qa,qa);
 % init damping force vector
Fem.system.fDamping = Fem.system.Damping * Fem.solver.sol.dx(qa);

% init global stiffness
Fem.system.Stiffness = sparse(Fem.triplets.i,Fem.triplets.j,           ...
    E(Fem.triplets.e).*Fem.triplets.k);  
Fem.system.Stiffness = Fem.system.Stiffness(qa,qa);

 % init global tangent stiffness
Fem.system.Tangent   = sparse(Fem.triplets.i,Fem.triplets.j,           ...
    E(Fem.triplets.e).*Fem.triplets.t); 
Fem.system.Tangent = Fem.system.Tangent(qa,qa);

   % init elastic force vector
Fem.system.fElastic  = sparse(Fem.triplets.i,1,E(Fem.triplets.e)       ...
    .*Fem.triplets.fi);  
Fem.system.fElastic = Fem.system.fElastic(qa);

 % init body force vector
Fem.system.fBody     = sparse(Fem.triplets.i,1,E(Fem.triplets.e)       ...
    .*Fem.triplets.fb);    
Fem.system.fBody = beta * Fem.system.fBody(qa);


    
end
