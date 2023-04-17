% 10-Apr-2023 08:10:35
% Auto-generated test script

% Initialize the test suite
% Add test cases here

clr;

sdf = sRectangle(0, 100, 0, 50);
msh = Mesh(sdf,'Quads',[30,5],'MaxIteration',150);
msh = msh.generate();

fem = Fem(msh,'TimeStep',1/10);
fem = fem.addMaterial(NeoHookean(1,0.4));
fem = fem.addSupport('left', [1, 1]);
%fem = fem.addSupport('right', [0,1]);

fem = fem.addDisplace('right', [150,0]);

fem = solveQuasiStaticFem(fem);


% qd  = fem.findNodes('right')*2-1;
% qa  = fem.system.Ia;

% C = eye(fem.Dim*fem.Mesh.NNode);
% C = C(qd,qa);

% Xd = 300;

% ndofq = numel(qa(qa));
% ndofp = numel(qd);

% N = 10;
% for ii = 1:N

%     beta = ii / N;
%     fem.options.loadingFactor = beta;

%     f = Inf;
%     jj = 1;

%     while norm(f) > 1e-3 && jj < 250
%         fem = fem.compute();
%         f = fem.system.fResidual + C.'*fem.solver.sol.u(qd);
%         K = fem.system.Tangent;
%         A = [K, C.'; C, zeros(ndofp)];
%         b = [f; fem.solver.sol.x(qd) - beta * Xd * ones(ndofp,1)];

%         dq = A\b;

%         fem.solver.sol.x(qa) = fem.solver.sol.x(qa) - dq(1:ndofq);
%         fem.solver.sol.u(qd) = fem.solver.sol.u(qd) - dq(ndofq+1:end);

%         log(ii,jj,f,round(mean(fem.solver.sol.x(qd))));

%         jj = jj + 1;
%     end

%     %fem.states.x(qa) = fem.solver.sol.x(qa);

%     h = msh;
%     h.Node = msh.Node + meshfield(fem, fem.solver.sol.x);
%     h.show();
%     drawnow;
% end

% axis on;

% function log(ii,jj,f,g)
%     if ii == 1 && jj == 1 
%         fprinttable({'itr', 'step','force residual','change'}, [ii, jj, norm(f),g],'open',true);
%     else
%         fprinttable({'itr', 'step', 'force residual','change'}, [ii, jj, norm(f),g], 'addrow',true,'open',true);
%     end
%     pause(.0);
% end
