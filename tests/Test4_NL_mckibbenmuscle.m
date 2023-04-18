% 10-Apr-2023 11:50:28
% Auto-generated test script

% Initialize the test suite
% Add test cases here

clf;
s = sRectangle(-5,5,0,40);
msh = Mesh(s,'Quads',[20,50]);
msh = msh.generate();

E = msh.findElements('Box',[-3,3,4,36]);
msh = msh.removeElements(E);

fem = Fem(msh,'TimeStep',1/10);
fem = fem.addMaterial(NeoHookean);
fem = fem.addMaterial(NeoHookean(100,0.3));
fem = fem.addMaterial(NeoHookean(10,0.3));

fem = fem.setMaterial('left',2);
fem = fem.setMaterial('right',2);
fem = fem.setMaterial(msh.findElements('box',[-5,5,0,4]),3);
fem = fem.setMaterial(msh.findElements('box',[-5,5,36,40]),3);

%showMaterialsFem(fem);

I = fem.findEdges('hole',2);
fem = fem.addPressure(I,150 * 1e-3);
fem = fem.addSupport('top',[1,1]);
qa = fem.system.Ia;

%fem.options.Display = @plt;

fem = solveQuasiStaticFem(fem);

% function plt(Fem)
%     h = Fem.Mesh;
%     h.Node = Fem.Mesh.Node + meshfield(Fem, Fem.solver.sol.x);
%     h.show();
%     axis([-50 50 -0 40]);
% end

% N = 50;
% for ii = 1:N

%     beta = ii / N;
%     fem.options.loadingFactor = beta;

%     f = Inf;
%     jj = 1;
%     U0 = fem.states.x(qa);

%     while norm(f) > 1e-3
%         fem = fem.compute();
%         f   = fem.system.fResidual;

%         dq  = fem.system.Tangent \ (fem.system.fResidual);

%         fem.solver.sol.x(qa) = fem.solver.sol.x(qa) - dq;

%         jj = jj + 1;
%     end

%     log(ii,jj,f);

%     fem.states.x(qa) = fem.solver.sol.x(qa);

%     h = msh;
%     h.Node = msh.Node + meshfield(fem, fem.solver.sol.x);
%     h.show();
%     drawnow;
% end

% showMaterialsFem(fem);

% function log(ii,jj,f)
%     if ii == 1 
%         fprinttable({'itr', 'k  ','force residual'}, [ii, jj, norm(f)],'open',true);
%     else
%         fprinttable({'itr', 'k  ', 'force residual'}, [ii, jj, norm(f)], 'addrow',true,'open',true);
%     end
% end