clr;

msh = Mesh('PneunetFine.png','BdBox',[0,120,0,20],'NElem',650);
msh = msh.generate();

msh.show();

fem = Fem(msh,'TimeStep',1/30);
fem = fem.addMaterial(NeoHookean(.2, 0.3));
fem = fem.addMaterial(NeoHookean(5, 0.33));

fem = fem.setMaterial(msh.findElements('box',[0,120,0,4]),2)

fem = fem.addPressure('allhole',@(t) 25 * 1e-3 * clamp(t,0,1));
fem = fem.addSupport('left',[1,1]);

showMaterialsFem(fem); pause(1);

fem.options.Display = @plt;

fem = solveDynamicFem(fem);

function plt(Fem)
    h = Fem.Mesh;
    h.Node = Fem.Mesh.Node + meshfield(Fem, Fem.solver.sol.x);
    h.show();
    axis([-50 120 -100 20]);
end
