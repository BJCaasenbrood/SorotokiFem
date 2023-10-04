function Fem = solve(Fem, varargin)
% SOLVE Solves a computational mechanics problem using the Finite Element Method.
%
%   Fem = solve(Fem, varargin) solves the given computational mechanics problem
%   using an adaptive Newton Rhapson method. The function takes in a Fem structure 
%   and optional input arguments specified by varargin.
%
%   Example:
%      Fem = Fem.solve()
%      Fem = Fem.solve('Display', @display);   % sets display function
%      Fem = Fem.solve('TimeStep', 1e-3);      % sets timesteps for solver
%
%   See also compute, simulate, optimize.

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
  
Ceq = [];
qc  = [];    
nc  = 0;
qa  = Fem.system.Ia;
nq  = numel(find(Fem.system.Ia));

% checks if the problem has displacement constraints
if isfield(Fem.system,'Displace')
    Fem = Fem.compute();
    qc  = Fem.system.Ic;
    nc  = numel(find(Fem.system.Ic));
    Ceq = Fem.system.cMatrix;
end

% checks if the problem is linear
if Fem.options.isNonlinear
    Fem.solver.Time  = 0;
else
    Fem.solver.Time = Fem.solver.TimeHorizon - 1e-12;
    Fem.solver.MaxIteration = 2;
end

NSteps = round(Fem.solver.TimeHorizon/Fem.solver.TimeStep);

% preallocate solutions
Fem.solver.sol.yout = zeros(NSteps,Fem.Dim * Fem.Mesh.NNode);
Fem.solver.sol.tout = zeros(NSteps,1);
Fem.solver.SubIteration = 1;

% Fem.log.info('Starting implicit Newton-Raphson method');
% Fem.log.hline();

if Fem.solver.isLog
    % progBar = Fem.log.progress(NSteps); 
    disp('|  Iter | Eval |   Residual   |  Time  |  Step  |');
end
    
% solves problem over [0, TimeHorizon]
while Fem.solver.Time < Fem.solver.TimeHorizon

    beta = Fem.solver.Time / Fem.solver.TimeHorizon + 1e-12;
    Fem.options.loadingFactor = beta;

    Fem.solver.Residual  = Inf;
    Fem.solver.Iteration = 1;

    x0 = Fem.solver.sol.x(qa);
    u0 = Fem.solver.sol.u(qc);

    % Newton-Raphson loop
    while norm(Fem.solver.Residual) > Fem.solver.RelTolerance && ...
        Fem.solver.Iteration < Fem.solver.MaxIteration 

        % if Fem.solver.Iteration < 5
        Fem = Fem.compute();

        % linear solve
        if isempty(Ceq)
            b = Fem.system.fResidual;
            A = Fem.system.Tangent;
        else
            b = [Fem.system.fResidual + Ceq.'*Fem.solver.sol.u(qc); ...
                 -Fem.system.cResidual];
         
            A = [Fem.system.Tangent, Ceq.'; Ceq, zeros(nc)];
        end

        % solve for the displacement increment
        % dA = decomposition(A,'qr');
        dfdq1 = A \ b;

        if Fem.solver.Iteration > 1
            minL = sqrt(1+th) * lam0;
            maxL = 0.5 * norm([x1;u1] - [x0;u0])/norm(dfdq1 - dfdq0);
            lam1 = clamp(min([minL, maxL]), 1e-6, Inf);
            th   = lam1/lam0;
        else
            lam0 = 1;
            lam1 = lam0;
            th   = +Inf;
        end

        x0 = Fem.solver.sol.x(qa);
        x1 = Fem.solver.sol.x(qa) - lam1 * dfdq1(1:nq);

        u0 = Fem.solver.sol.u(qc);
        u1 = Fem.solver.sol.u(qc) - lam1 * dfdq1(nq+(1:nc));

        Fem.solver.Residual  = b;
        Fem.solver.Iteration = Fem.solver.Iteration + 1;
        Fem.solver.sol.u(qc) = u1;
        Fem.solver.sol.x(qa) = x1;
        lam0  = lam1;
        dfdq0 = dfdq1;

        if Fem.solver.isLog
            fprintf(repmat('\b',1,80));
            s=sprintf(' %5.0f %5.0f %14.2e %9.2g %8.2g', ...
                Fem.solver.SubIteration, Fem.solver.Iteration-1, norm(Fem.solver.Residual), Fem.solver.Time,lam1); 
            fprintf(s);
        end
    end

    step = Fem.solver.SubIteration;
    Fem.solver.sol.yout(step,:) = Fem.solver.sol.x;
    Fem.solver.sol.tout(step)   = Fem.solver.Time;
    Fem.solver.SubIteration = Fem.solver.SubIteration + 1;

    Fem.solver.Time = clamp(Fem.solver.Time + Fem.solver.TimeStep,...
        0,Fem.solver.TimeHorizon);

    if Fem.solver.isLog && Fem.solver.Time < Fem.solver.TimeHorizon
        % progBar(step);
    end

    if ~isempty(Fem.options.Display)
        Fem.options.Display(Fem);
        drawnow;
    end
end

% fprintf('\n');
% Fem.log.hline();
% Fem.log.info('Generalized Newton-Raphson completed');
% Fem.log.info('States written to:');
% Fem.log.list('',{'states = *.solver.sol.yout'});

end