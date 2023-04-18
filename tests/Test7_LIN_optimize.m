% 10-Apr-2023 08:10:35
% Auto-generated test script

% Initialize the test suite
% Add test cases here
clr;

sdf = sRectangle(0, 5, 0, 2);
msh = Mesh(sdf,'Quads',[80,40],'MaxIteration',150);
msh = msh.generate();

fem = Fem(msh,'TimeStep',0.3,'SpatialFilterRadius',0.1,...
          'Interpolation','SIMP-H','MaxChange',Inf,'Penal',1);

fem = fem.addMaterial(NeoHookean);
fem = fem.addSupport('left', [1, 0]);
fem = fem.addSupport('se',[0, 1]);
fem = fem.addLoad(fem.findNodes('Location',[0,2],1),[0,1e-3]);
fem = fem.addLoad(fem.findNodes('Location',[0,0],1),[-1e-3,0]);

sdf = sCircle(0.5,[2.5,1]);
I   = sdf.intersect(msh.Center);

fem.topology.sol.x(I) = 1e-6;
fem.options.isNonlinear = false;
fem.options.LineStyle = 'none';
fem.options.Display = @plt;

fem = solveOptimizationFem(fem);

function plt(Fem)
    cla;
    showInfillFem(Fem);
end
