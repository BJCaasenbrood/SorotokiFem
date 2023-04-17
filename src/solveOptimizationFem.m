function Fem = solveQuasiStaticFem(Fem, varargin)

Fem.solver.Time = 0;
qa = Fem.system.Ia;

while Fem.topology.Iteration < Fem.topology.MaxIteration
        
        % update material field
        [~,dEdy,~,dVdy] = MaterialField(Fem);
         
        % set load factor
        % Fem.OptFactor = sigmoid(5*Fem.IterationMMA/Fem.MaxIterationMMA,...
        %     Fem.SigmoidFactor);
       
        % solve nonlinear finite elements
        Fem = Fem.solve();
        
        % evaluate objective and constraints
        [f,dfdE,dfdV] = ObjectiveFunction(Fem);
        [g,dgdE,dgdV] = ConstraintFunction(Fem);
    
        if ~isempty(Fem.Source)
            dfdE(Fem.Source(:,1)) = -1e-4*Fem.ChangeMax;
        end
        
        % compute design sensitivities
        dfdz = Fem.SpatialFilter'*(dEdy.*dfdE + dVdy.*dfdV);
        dgdz = Fem.SpatialFilter'*(dEdy.*dgdE + dVdy.*dgdV);
        
        %compute design variable
        [Fem,ZNew] = UpdateSchemeMMA(Fem,f,dfdz,g,dgdz);
        
        % determine material change
        Fem.Change  = clamp(ZNew - Fem.Density,-Fem.ChangeMax,Fem.ChangeMax);

    % beta = Fem.solver.Time / Fem.solver.TimeHorizon;
    % Fem.options.loadingFactor = beta;

    % Fem.solver.Residual  = Inf;
    % Fem.solver.Iteration = 1;

    % x0 = Fem.solver.sol.x(qa);

    % while norm(Fem.solver.Residual) > Fem.solver.RelTolerance && ...
    %     Fem.solver.Iteration < Fem.solver.MaxIteration 

    %     Fem = Fem.compute();
    %     f   = Fem.system.fResidual;

    %     % linear solve
    %     dfdq1 = Fem.system.Tangent \ (Fem.system.fResidual);

    %     if Fem.solver.Iteration > 1
    %         minL = sqrt(1+th) * lam0;
    %         maxL = 0.5 * norm(x1 - x0)/norm(dfdq1 - dfdq0);
    %         lam1 = clamp(min([minL, maxL]),1,Inf);
    %         th   = lam1/lam0;
    %     else
    %         lam0 = 1;
    %         lam1 = lam0;
    %         th   = +Inf;
    %     end

    %     Fem.solver.Residual  = f;
    %     Fem.solver.Iteration = Fem.solver.Iteration + 1;

    %     dfdq0 = dfdq1;
    %     x0 = Fem.solver.sol.x(qa);
    %     x1 = Fem.solver.sol.x(qa) - lam1 * dfdq1;

    %     Fem.solver.sol.x(qa) = x1;
    %     lam0 = lam1;

        log(Fem.solver.Time + 1e-6, ...
            Fem.solver.Iteration,...
            norm(f),...
            lam1 + 1e-6);
    end


    % Fem.solver.Time = clamp(Fem.solver.Time + Fem.solver.TimeStep,...
    %     0,Fem.solver.TimeHorizon);
    % Fem.states.x(qa) = Fem.solver.sol.x(qa);
end

end

function log(ii,jj,f,g)
    if ii < 1e-3 
        fprinttable({'time', 'step','force residual','lambda'}, [ii, jj, f,g],'open',true);
    else
        fprinttable({'time', 'step', 'force residual','lambda'}, [ii, jj, f,g], 'addrow',true,'open',true);
    end
    pause(.0);
end


% %---------------------------------------------- topology optimization solve
% function Fem = optimize(Fem)
 
%     if Fem.ShowProcess
%         showInformation(Fem,'TopologyOptimization');
%     end    
        
%     Fem.SolverPlot     = false;
%     Fem.IterationMMA   = 0;
%     Fem.SolverStartMMA = true;
%     flag               = true;
%     Visual             = 'ISO';
%     PenalUpdate        = round(Fem.MaxIterationMMA/5);
%     Fem.SpatialFilter  = GenerateRadialFilter(Fem,Fem.FilterRadius);
    
%     % draw initial visual
%     Fem.show(Visual); drawnow;
    
%     while flag
    
%         % update increment
%         Fem.IterationMMA = Fem.IterationMMA + 1;
        
%         % update material field
%         [~,dEdy,~,dVdy] = MaterialField(Fem);
         
%         % set load factor
%         Fem.OptFactor = sigmoid(5*Fem.IterationMMA/Fem.MaxIterationMMA,...
%             Fem.SigmoidFactor);
       
%         % solve nonlinear finite elements
%         Fem = Fem.solve();
        
%         % evaluate objective and constraints
%         [f,dfdE,dfdV] = ObjectiveFunction(Fem);
%         [g,dgdE,dgdV] = ConstraintFunction(Fem);
    
%         if ~isempty(Fem.Source)
%             dfdE(Fem.Source(:,1)) = -1e-4*Fem.ChangeMax;
%         end
        
%         % compute design sensitivities
%         dfdz = Fem.SpatialFilter'*(dEdy.*dfdE + dVdy.*dfdV);
%         dgdz = Fem.SpatialFilter'*(dEdy.*dgdE + dVdy.*dgdV);
        
%         %compute design variable
%         [Fem,ZNew] = UpdateSchemeMMA(Fem,f,dfdz,g,dgdz);
        
%         % determine material change
%         Fem.Change  = clamp(ZNew - Fem.Density,-Fem.ChangeMax,Fem.ChangeMax);
    
%         Fem.Density = Fem.Density + Fem.Change;
        
%         % evaluate fitness
%         Fem.Objective  = vappend(Fem.Objective,f);       
%         Fem.Constraint = vappend(Fem.Constraint,g);
        
%         % check convergence
%         [flag,Fem] = CheckConvergenceOpt(Fem);
        
%         % save data
%         Fem = FiniteElementDataLog(Fem);
        
%         % draw visual
%         if (mod(Fem.IterationMMA,1) == 0) || Fem.Nonlinear
%             Fem.show(Visual); drawnow;
%             colormap(turbo);
%         end
        
%         if mod(Fem.IterationMMA,PenalUpdate) == 0
%            Fem.Penal = clamp(Fem.Penal + 1,1,5);
%         end
%     end
    
%     Fem.SolverStartMMA = false;
    
%     end