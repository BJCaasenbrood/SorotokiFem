clr;

sdf = sRectangle(0, 100, 0, 3);
msh = Mesh(sdf,'Quads',[50,2],'MaxIteration',150);
msh = msh.generate();

fem = Fem(msh,'TimeStep',1/750);
fem = fem.addMaterial(NeoHookean(1,0.49));
fem = fem.addSupport('left', [1, 1]);
fem = fem.addContact(sCircle(10,[70,-25]));
fem = fem.addGravity();

fem.options.Display = @plt;

fem = solveDynamicFem(fem);


function plt(Fem)
    h = Fem.Mesh;
    h.Node = Fem.Mesh.Node + meshfield(Fem, Fem.solver.sol.x);
    h.show();
    axis([-50 120 -100 20]);
end