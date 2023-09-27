function Fem = simulate(Fem, varargin)
% SOLVE simulates dynamic computational mechanics problem using the Finite Element Method.
%
%   Fem = simulate(Fem, varargin) solves the dynamical computational mechanics problem
%   using an adaptive Beta-Newmark method. The function takes in a Fem structure 
%   and optional input arguments specified by varargin.
%
%   Example:
%      Fem = Fem.simulate()
%      Fem = Fem.simulate('Display', @display);   % sets display function
%      Fem = Fem.simulate('TimeStep', 1e-3);      % sets timesteps for solver
%
%   See also compute, simulate, optimize, update.

for ii = 1:2:length(varargin)
    if isprop(Fem.options,varargin{ii})
        Fem.options.(varargin{ii}) = varargin{ii+1};
    elseif isprop(Fem.solver,varargin{ii})
        Fem.solver.(varargin{ii}) = varargin{ii+1};
    elseif isprop(Fem.topology,varargin{ii})
        Fem.topology.(varargin{ii}) = varargin{ii+1};
    else
        Fem.(varargin{ii}) = varargin{ii+1};
    end
end

rho = 0;  % must betwween [0,1]
alphaM = (2*rho - 1)/(rho + 1);
alphaF = (rho)/(rho + 1);

gamma = 0.5 - alphaM + alphaF;
beta  = (0.25)*(1 - alphaM + alphaF)^2;

Fem.solver.Time = 0;
if ~isfield(Fem.system,'Ia')
    Fem.system.Ia = 1:Fem.Dim * Fem.Mesh.NNode;
end

if Fem.solver.isLog
    NSteps = round(Fem.solver.TimeHorizon/Fem.solver.TimeStep);
    % progBar = ProgressBar(NSteps,'Title', ' ');
end

qa = Fem.system.Ia; 

if ~isfield(Fem.solver.sol,'ddx')

    Fem.solver.sol.ddx = Fem.solver.sol.dx * 0;
    Fem = Fem.compute();

    A = Fem.system.Mass;
    b = -Fem.system.fResidual - Fem.system.Damping * ...
        Fem.solver.sol.dx(qa);

    Fem.solver.ddx(qa) = A \ b;
end

while Fem.solver.Time < Fem.solver.TimeHorizon

    dt  = Fem.solver.TimeStep;
    Fem.solver.Residual  = Inf;
    Fem.solver.Iteration = 1;

    Fem.options.loadingFactor = smoothstep(Fem.solver.Time * 5); 

    x0   = Fem.solver.sol.x(qa);
    dx0  = Fem.solver.sol.dx(qa);
    ddx0 = Fem.solver.sol.ddx(qa);
    ddxf = Fem.solver.sol.ddx(qa);
    tf   = Fem.solver.Time;

    while norm(Fem.solver.Residual) > Fem.solver.RelTolerance && ...
        Fem.solver.Iteration < Fem.solver.MaxIteration 

        xf  = x0  + dt * dx0 + dt*dt * ((.5-beta)*ddxf + beta*ddx0);
        dxf = dx0 + dt * ((1-gamma)*ddxf + gamma*ddx0);

        x1   = (1-alphaF)*xf + alphaF*x0;
        dx1  = (1-alphaF)*dxf + alphaF*dx0;     
        ddx1 = (1-alphaM)*ddx0 + alphaM*ddxf;      

        Fem.solver.sol.x(qa)   = x1;
        Fem.solver.sol.dx(qa)  = dx1;
        Fem.solver.sol.ddx(qa) = ddx1;
        Fem.solver.Time = tf + dt - alphaF*dt;

        Fem = Fem.compute();

        A = (1-alphaF) * beta * dt*dt * Fem.system.Tangent + ...
            (1-alphaF) * gamma * dt * Fem.system.Damping + ...
            (1-alphaM) * Fem.system.Mass;

        b = Fem.system.Mass * ddx1 + Fem.system.fResidual;
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


        if Fem.solver.isLog
            if Fem.solver.Iteration == 2,
                disp(' -----------------------------------------');
                disp(' |  Iter | Eval |   f(x)   |  Step-size  |');
                disp(' -----------------------------------------');
            end
            s=sprintf(' %5.0f  %5.0f     %5.3e %13.5g    ', ...
                Fem.solver.SubIteration, Fem.solver.Iteration-1, norm(Fem.solver.Residual),lam1); 
            disp(s);
        end

    end

    % if Fem.solver.isLog
    %     progBar([], [], []);
    % end

    Fem.solver.Time = clamp(tf + Fem.solver.TimeStep,...
        0,Fem.solver.TimeHorizon);

    Fem.solver.sol.dx(qa) = dx0 + dt * ((1-gamma)*ddxf + gamma * ddx0);    
    Fem.solver.sol.x(qa)  = x0 + dt * dx0 + dt*dt * ((0.5 - beta)*ddxf + beta * ddx0);

    if ~isempty(Fem.options.Display)
        Fem.options.Display(Fem);
        drawnow;
    end
end

% if Fem.solver.isLog
%     progBar.release();
%     fprintf('\n');
% end

end
