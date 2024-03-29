function Fem = update(Fem, varargin)

if nargin > 1
    Fem.solver.TimeStep = varargin{1};
end    

rho    = 0.25;
alphaM = (2*rho - 1)/(rho + 1);
alphaF = (rho)/(rho + 1);

gamma = 0.5 - alphaM + alphaF;
beta  = 0.25*(1 - alphaM + alphaF)^2;

if ~isfield(Fem.system,'Ia')
    Fem.system.Ia = 1:Fem.Dim * Fem.Mesh.NNode;
end

qa = Fem.system.Ia;

if Fem.solver.Time < Fem.solver.TimeStep  
    if ~isfield(Fem.solver.sol,'ddx')
        
        Fem.solver.sol.ddx = Fem.solver.sol.dx * 0;
        Fem = Fem.compute();
        
        A = Fem.system.Mass;
        b = -Fem.system.fResidual - Fem.system.Damping * ...
            Fem.solver.sol.dx(qa);
        
        Fem.solver.ddx(qa) = A \ b;
    end

    Fem.options.loadingFactor = 1;
end

dt  = Fem.solver.TimeStep;
Fem.solver.Residual  = Inf;
Fem.solver.Iteration = 1;

x0   = Fem.solver.sol.x(qa);
dx0  = Fem.solver.sol.dx(qa);
ddx0 = Fem.solver.sol.ddx(qa);
ddxf = Fem.solver.sol.ddx(qa);
tf   = Fem.solver.Time;

while norm(Fem.solver.Residual) > Fem.solver.RelTolerance && ...
        Fem.solver.Iteration < Fem.solver.MaxIteration
    
    dxf = dx0 + dt * ((1-gamma)*ddxf + gamma*ddx0);
    xf  = x0  + dt * dx0 + 0.5 * dt*dt * ...
        ((1-2*beta)*ddxf + 2*beta*ddx0);
    
    x1   = (1-alphaF)*xf + alphaF*x0;
    dx1  = (1-alphaF)*dxf + alphaF*dx0;
    ddx1 = (1-alphaF)*ddx0 + alphaF*ddxf;
    
    Fem.solver.sol.x(qa)   = x1;
    Fem.solver.sol.dx(qa)  = dx1;
    Fem.solver.sol.ddx(qa) = ddx1;
    Fem.solver.Time        = tf + dt - alphaF*dt;
    
    Fem = Fem.compute();
    
    A = (1-alphaF) * beta * dt*dt * Fem.system.Tangent + ...
        (1-alphaF) * gamma * dt * Fem.system.Damping + ...
        (1-alphaM) * Fem.system.Mass;
    
    b = Fem.system.Mass * ddx1 + Fem.system.fResidual;
    
    % linear solve
    dfdq1 = A \ b;
    
    if Fem.solver.Iteration > 1
        minL = sqrt(1+th) * lam0;
        maxL = norm(ddxf_ - ddx0)/norm(dfdq1 - dfdq0);
        lam1 = clamp(min([minL, maxL]),0,Inf);
    else
        lam0 = 1;
        lam1 = lam0;
        th   = +Inf;
    end
    
    Fem.solver.Residual  = b;
    Fem.solver.Iteration = Fem.solver.Iteration + 1;
    
    ddxf_ = ddx0;
    ddx0  = ddx0 - lam1 * dfdq1;
    dfdq0 = dfdq1;
end

% Fem.solver.Time = clamp(tf + Fem.solver.TimeStep,...
%     0,Fem.solver.TimeHorizon);

Fem.solver.sol.dx(qa) = dx0 + dt * ((1-gamma)*ddxf + gamma * ddx0);
Fem.solver.sol.x(qa)  = x0 + dt * dx0 + 0.5 * dt*dt * ((1 - 2*beta)*ddxf + 2 * beta * ddx0);

end
