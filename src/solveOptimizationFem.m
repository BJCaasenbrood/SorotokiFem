function Fem = solveOptimizationFem(Fem, varargin)

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

while Fem.topology.Iteration < Fem.topology.MaxIteration
        
        % update material field
        [~,dEdy,~,dVdy] = materialFieldFem(Fem);
         
        % set load factor
        % Fem.OptFactor = sigmoid(5*Fem.IterationMMA/Fem.MaxIterationMMA,...
        %     Fem.SigmoidFactor);
       
        % solve nonlinear finite elements
        Fem = solveQuasiStaticFem(Fem);
        
        % evaluate objective and constraints
        [f,dfdE,dfdV] = assembleTopoObjective(Fem);
        [g,dgdE,dgdV] = assembleTopoConstraint(Fem);
      
        % compute design sensitivities
        Filter = Fem.topology.SpatialFilter;
        dfdz = Filter'*(dEdy.*dfdE + dVdy.*dfdV);
        dgdz = Filter'*(dEdy.*dgdE + dVdy.*dgdV);
        
        if Fem.topology.Iteration > 1, x1 = x1; else,  x1 = 0; end
        if Fem.topology.Iteration > 2, x2 = x2; else,  x2 = 0; end
        
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Fem,zNew] = UpdateSchemeMMA(Fem,f,dfdz,g,dgdz)  

iter = Fem.topology.Iteration; 
N    = length(Fem.topology.sol.x);
M    = 1;

if isempty(Fem.fnorm)
    Fem.fnorm = abs(norm(dfdz)); 
end

f0val  = f;
df0dx  = dfdz;
df0dx2 = 0*df0dx;
fval   = g;
dfdx   = dgdz;
dfdx2  = dgdz*0;

xval   = Fem.topology.sol.x;
xmin   = zeros(N,1);
xmax   = ones(N,1);
A0 = 1; 
A  = 0; 
C  = 10000*ones(M,1); 
D = 0;

if iter == 1, Fem.low = xmin; Fem.upp = xmax; end
if iter > 1, Fem.xold1 = Fem.xold1; else, Fem.xold1 = 0; end
if iter > 2, Fem.xold2 = Fem.xold2; else, Fem.xold2 = 0; end

[xmma,~,~,~,~,~,~,~,~,Fem.low,Fem.upp] = mmasub(M,N,iter,xval,xmin,...
    xmax,Fem.xold1,Fem.xold2,f0val,df0dx,df0dx2,fval,dfdx,dfdx2,Fem.low,...
    Fem.upp,A0,A,C,D);

dx   = clamp(xmma - xval,-Fem.ChangeMax*15,Fem.ChangeMax*15);            
zNew = xval + dx;

if iter >= 1, 
    Fem.xold1 = xval; 
end

if iter >= 2, 
    Fem.xold2 = Fem.xold1;
 end

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