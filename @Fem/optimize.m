function Fem = optimize(Fem, varargin)

% ensure that opt-log is not over run fem-log
Fem.solver.Display = false;   
Fem.solver.isLog   = false;   

Hmax = Fem.topology.MaxChange;
xmin = zeros(Fem.Mesh.NElem,1);
xmax = ones(Fem.Mesh.NElem,1);
upper = xmax;
lower = xmin;
A0   = 1; 
A    = 0; 
C    = 10000; 
D    = 0;

%Fem.solver.Time = 0;
%qa = Fem.system.Ia;

assert(Fem.materials.NMat == 1, 'Currently, single material optimization is only supported.');

x0 = Fem.topology.sol.x;
x1 = Fem.topology.sol.x;
x2 = Fem.topology.sol.x;

progBar = ProgressBar(Fem.topology.MaxIteration,'Title', ' ');

while Fem.topology.Iteration <= Fem.topology.MaxIteration
        
        % update material field
        [~,dEdy,~,dVdy] = materialFieldFem(Fem);
         
        % set load factor
        % Fem.OptFactor = sigmoid(5*Fem.IterationMMA/Fem.MaxIterationMMA,...
        %     Fem.SigmoidFactor);
       
        % solve nonlinear finite elements
        Fem = solve(Fem);
        
        % evaluate objective and constraints
        [f,dfdE,dfdV] = assembleTopoObjective(Fem);
        [g,dgdE,dgdV] = assembleTopoConstraint(Fem);
      
        % compute design sensitivities
        Filter = Fem.topology.SpatialFilter;
        dfdz = Filter'*(dEdy.*dfdE + dVdy.*dfdV);
        dgdz = Filter'*(dEdy.*dgdE + dVdy.*dgdV);
        
        if Fem.topology.Iteration <= 1, x1 = 0; end
        if Fem.topology.Iteration <= 2, x2 = 0; end
        
        [xmma,~,~,~,~,~,~,~,~,lower,upper] = ...
            mmasub(1,Fem.Mesh.NElem,Fem.topology.Iteration,...
            x0,xmin,xmax,x1,x2,...
            f,dfdz,0*dfdz,g,dgdz,0*dgdz,...
            lower,upper,A0,A,C,D);

        if Fem.topology.Iteration >= 1, x1 = x0; end
        if Fem.topology.Iteration >= 2, x2 = x1; end
        x0 = x0 + clamp(xmma-x0,-Hmax,Hmax);

        Fem.topology.sol.x     = x0;
        Fem.topology.Iteration = Fem.topology.Iteration + 1;
        Fem.topology.Penal = clamp(Fem.topology.Penal +  ...
            Fem.topology.PenalStep, 0, Fem.topology.MaxPenal);

        if ~isempty(Fem.options.Display)
            Fem.options.Display(Fem);
            drawnow;
        end

        progBar([], [], []);
end

progBar.release();
fprintf('\n');

end
