function Fem = solveQuasiStaticFem(Fem, varargin)
  
Ceq = [];
qc  = [];    
nc  = 0;
qa  = Fem.system.Ia;
nq  = numel(find(Fem.system.Ia));

if isfield(Fem.system,'Displace')
    Fem = Fem.compute();
    qc  = Fem.system.Ic;
    nc  = numel(find(Fem.system.Ic));
    Ceq = Fem.system.cMatrix;
end

if Fem.options.isNonlinear
    Fem.solver.Time  = 0;
else
    Fem.solver.Time = Fem.solver.TimeHorizon - 1e-12;
    Fem.solver.MaxIteration = 2;
end

if Fem.solver.isLog
    NSteps = round(Fem.solver.TimeHorizon/Fem.solver.TimeStep);
    progBar = ProgressBar(NSteps,'Title', 'Solve static FEM');
end

    % preallocate solutions
    Fem.solver.sol.yout = zeros(NSteps,Fem.Dim * Fem.Mesh.NNode);
    Fem.solver.sol.tout = zeros(NSteps,1);
    Fem.solver.SubIteration = 1;
    
while Fem.solver.Time < Fem.solver.TimeHorizon

    beta = Fem.solver.Time / Fem.solver.TimeHorizon;
    Fem.options.loadingFactor = beta;

    Fem.solver.Residual  = Inf;
    Fem.solver.Iteration = 1;

    x0 = Fem.solver.sol.x(qa);
    u0 = Fem.solver.sol.u(qc);

    while norm(Fem.solver.Residual) > Fem.solver.RelTolerance && ...
        Fem.solver.Iteration < Fem.solver.MaxIteration 

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

        dfdq1 = A \ b;

        if Fem.solver.Iteration > 1
            minL = sqrt(1+th) * lam0;
            maxL = 0.5 * norm([x1;u1] - [x0;u0])/norm(dfdq1 - dfdq0);
            lam1 = clamp(min([minL, maxL]),1,1);
            th   = lam1/lam0;
        else
            lam0 = 1;
            lam1 = lam0;
            th   = +Inf;
        end

        % if  Fem.solver.Display && (mod(Fem.solver.Iteration,Fem.solver.DisplayAtEvery) == 0 || ...
        %     Fem.solver.Iteration == 1)

        %     log(Fem.solver.Time + 1e-6, Fem.solver.Iteration,...
        %         norm(b), lam1 + 1e-6);
        % end

        Fem.solver.Residual  = b;
        Fem.solver.Iteration = Fem.solver.Iteration + 1;

        dfdq0 = dfdq1;
        x0 = Fem.solver.sol.x(qa);
        x1 = Fem.solver.sol.x(qa) - lam1 * dfdq1(1:nq);

        u0 = Fem.solver.sol.u(qc);
        u1 = Fem.solver.sol.u(qc) - lam1 * dfdq1(nq+(1:nc));

        Fem.solver.sol.u(qc) = u1;
        Fem.solver.sol.x(qa) = x1;
        lam0 = lam1;

    end

    if Fem.solver.isLog
        progBar([], [], []);
    end

    step = Fem.solver.SubIteration;
    Fem.solver.sol.yout(step,:) = Fem.solver.sol.x;
    Fem.solver.sol.tout(step)   = Fem.solver.Time;
    Fem.solver.SubIteration = Fem.solver.SubIteration + 1;

    Fem.solver.Time = clamp(Fem.solver.Time + Fem.solver.TimeStep,...
        0,Fem.solver.TimeHorizon);

    if ~isempty(Fem.options.Display)
        Fem.options.Display(Fem);
        drawnow;
    end
end

if Fem.solver.isLog
    fprintf('\n');
end

end

% function log(ii,jj,f,g)
%     if ii < 1e-3 
%         fprinttable({'time', 'step','force residual','lambda'}, [ii, jj, f,g],'open',true);
%     else
%         fprinttable({'time', 'step', 'force residual','lambda'}, [ii, jj, f,g], 'addrow',true,'open',true);
%     end
%     pause(.0);
% end

