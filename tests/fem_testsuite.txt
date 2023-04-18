clr;
% f1 = @(x) sqrt(x(:,1).^2 + x(:,2).^2) - 1.0;
%
% B = [-2,2,-2,2];
% sdf1 = Sdf(f1,'BdBox',B);
sdf = sRectangle(0, 100, 0, 20);
msh = Mesh(sdf,'NElem', 120,'MaxIteration',50);
% msh = Mesh('PneunetFine.png','BdBox',[0,120,0,20]);
% msh = msh.set('SimplificationTolerance',1);
msh = msh.generate();

% %%
fem = Fem(msh);
fem = fem.addMaterial(NeoHookean);

fem = fem.addSupport('left', [1, 1]);
fem = fem.addTendon('right', [0, 1e-2]);
%fem.system.Gravity = [0;-9810];

qa = fem.system.Ia;

N = 120;

for ii = 1:N

    beta = ii / N;
    fem.options.loadingFactor = beta;

    f = Inf;
    jj = 1;
    %U0 = fem.states.x(qa);

    while norm(f) > 1e-3
        fem = fem.compute();
        f   = fem.system.fResidual;
        dq  = fem.system.Tangent \ (fem.system.fResidual);

        %fprintf('i=%f, r=%f\n', ii, norm(f));
        %T = evalc('fprinttable([{'itr','residual'}], [ii, norm(f)])');

        % [U,V] = eigs(fem.system.global.Mass,...
        %              fem.system.global.Tangent,5);

        %x = fem.solver.sol.x;
        fem.solver.sol.x(qa) = fem.solver.sol.x(qa) - dq;
        log(ii,jj,f);

        jj = jj + 1;
    end

    fem.states.x(qa) = fem.solver.sol.x(qa);

    h = msh;
    h.Node = msh.Node + meshfield(fem, fem.solver.sol.x);
    h.show('Un');
    drawnow;
end

fem.show('Un');


function log(ii,jj,f)
    if ii == 1 && jj == 1 
        fprinttable({'itr', 'k','force residual'}, [ii, jj, norm(f)],'open',true);
    else
        fprinttable({'itr', 'k', 'force residual'}, [ii, jj, norm(f)], 'addrow',true,'open',true);
    end
end
