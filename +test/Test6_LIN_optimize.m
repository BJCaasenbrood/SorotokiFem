% 10-Apr-2023 08:10:35
% Auto-generated test script

% Initialize the test suite
% Add test cases here
clr;

sdf = sRectangle(0, 10, 0, 2);
msh = Mesh(sdf,'Quads',[100,20],'MaxIteration',50);
msh = msh.generate();

fem = Fem(msh,'TimeStep',0.3,'SpatialFilterRadius',0.2,...
          'Interpolation','SIMP-H','MaxChange',Inf,'Penal',3);

fem = fem.addMaterial(NeoHookean);
fem = fem.addSupport('se',[1, 1]);
fem = fem.addSupport('sw',[1, 1]);
fem = fem.addLoad('bottommid',[0,1e-1]);

%fem.topology.sol.x(I) = 1e-6;
fem.options.isNonlinear = false;
fem.options.LineStyle = 'none';
fem.options.Display = @plt;

fem = solveOptimizationFem(fem);

function plt(Fem)
    cla;
    showInfillFem(Fem);
end
