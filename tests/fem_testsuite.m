clr;
% f1 = @(x) sqrt(x(:,1).^2 + x(:,2).^2) - 1.0;
%
% B = [-2,2,-2,2];
% sdf1 = Sdf(f1,'BdBox',B);
sdf = sRectangle(0, 100, 0, 5);
msh = Mesh(sdf, 'NElem', 1e3, 'MaxIteration', 50);
% msh = Mesh('PneunetFine.png','BdBox',[0,120,0,20]);
% msh = msh.set('SimplificationTolerance',1);
msh = msh.generate();

msh.show();
drawnow;

% %%
fem = Fem(msh);
fem = fem.addSupport('left', [1, 1]);
%fem = fem.addLoad('right', [0, 1e-3]);
fem.system.Gravity = [0;-9810];

% x = fem.solver.sol.x;

% while norm(gradobj(x0))>0.001
%
%     X = [X; x0];
%     linedir = -gradobj(x0);
%
%     fline = @(a) fobj(x0+a*linedir);
%     a =fminsearch(fline,0);
%
%     x0 = x0 + a*linedir;
%     k = k+1;
% end
% options = optimoptions('fminunc','Algorithm','trust-region',...
%     'SpecifyObjectiveGradient',true,'Display','iter');
%
% problem.options = opts;
% problem.x0 = fem.solver.sol.x;
% problem.objective = @(x) femsys(x,fem);
% problem.solver = 'fminunc';
%
% x = fminunc(problem);

% function [f, g, H] = femsys(x,Fem)
% % Calculate objective f
% Fem = Fem.compute(x);
% f = Fem.system.Potential;
% %
% % if nargout > 1 % gradient required
%      g = Fem.system.fElastic;
% %
% %     if nargout > 2 % Hessian required
% %         H = fem.system.Tangent;
% %     end
% %
% % end
% end1

N = 5;
qa = fem.system.Ia;

for ii = 1:N

    fem.options.loadingFactor = ii / N;

    f = Inf;
    jj = 1;
    U0 = fem.states.x(qa);

    while norm(f) > 0.1
        fem = fem.compute();
        f  = fem.system.fResidual;
        dq = fem.system.Tangent \ fem.system.fResidual;

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
        fprinttable([{'itr', 'k','force residual'}], [ii, jj, norm(f)],'open',true);
    else
        fprinttable([{'itr', 'k', 'force residual'}], [ii, jj, norm(f)], 'addrow',true,'open',true);
    end
end
