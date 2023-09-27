% 10-Apr-2023 08:10:35
% Auto-generated test script

% Initialize the test suite
% Add test cases here
clr;

msh = Mesh(sRectangle(10, 2),'Quads',[100,20]);
msh = msh.generate();

fem = Fem(msh,'SpatialFilterRadius',0.3);

fem = fem.addMaterial(NeoHookean);
fem = fem.addSupport('se',[1, 1]);
fem = fem.addSupport('sw',[1, 1]);
fem = fem.addLoad('bottommid',[0,1]);

fem.options.isNonlinear = false;
fem.options.LineStyle = 'none';
fem.options.Display = @plt;

fem = fem.optimize;

function plt(Fem)
    cla;
    showInfillFem(Fem);
end


