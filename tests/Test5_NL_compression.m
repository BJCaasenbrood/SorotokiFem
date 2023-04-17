% 10-Apr-2023 08:10:35
% Auto-generated test script

% Initialize the test suite
% Add test cases here

clr;

sdf = sRectangle(0, 10, 0, 10);
msh = Mesh(sdf,'Quads',[10,10],'MaxIteration',150);
%msh = Mesh(sdf,'NElem',100,'MaxIteration',150);
msh = msh.generate();

msh = msh.removeElements(1:2);
msh.show();
pause();

fem = Fem(msh);
fem = fem.addMaterial(NeoHookean(1,0.49));
fem = fem.addSupport('top', [1, 0]);
fem = fem.addSupport('bottom', [1,1]);
qd  = fem.findNodes('top')*2;
qa  = fem.system.Ia;

C = eye(fem.Dim*fem.Mesh.NNode);
C = C(qd,qa);

Xd = -4;

ndofq = numel(qa(qa));
ndofp = numel(qd);

N = 100;
for ii = 1:N

    beta = ii / N;
    fem.options.loadingFactor = beta;

    f = Inf;
    jj = 1;

    while norm(f) > 1e-10 && jj < 50
        fem = fem.compute();
        f = fem.system.fResidual + C.'*fem.solver.sol.u(qd);
        K = fem.system.Tangent;
        A = [K, C.'; C, zeros(ndofp)];
        b = [f; fem.solver.sol.x(qd) - beta * Xd * ones(ndofp,1)];

        dq = A\b;

        fem.solver.sol.x(qa) = fem.solver.sol.x(qa) - dq(1:ndofq);
        fem.solver.sol.u(qd) = fem.solver.sol.u(qd) - dq(ndofq+1:end);

        jj = jj + 1;
    end

    log(ii,jj,f,round(mean(fem.solver.sol.x(qd))));

    fem.states.x(qa) = fem.solver.sol.x(qa);

    h = msh;
    h.Node = msh.Node + meshfield(fem, fem.solver.sol.x);
    h.show();
    axis([-5,15,0,10]);
    drawnow;
end

axis on;

function log(ii,jj,f,g)
    if ii == 1 
        fprinttable({'itr', 'step','force residual','change'}, [ii, jj, norm(f),g],'open',true);
    else
        fprinttable({'itr', 'step', 'force residual','change'}, [ii, jj, norm(f),g], 'addrow',true,'open',true);
    end
    pause(.0);
end
