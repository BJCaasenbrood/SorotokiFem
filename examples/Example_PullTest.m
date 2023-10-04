clr;
sdf = sRectangle(0,1,0,1);
msh = Mesh(sdf,'Quads',[12,4]);
msh = msh.generate();

fem = Fem(msh,'TimeStep',1/30);
fem = fem.addSupport('Left',[1,1]);
fem = fem.addMaterial(NeoHookean(0.1,0.49));

fem = fem.addDisplace('right',[5,0]);

fem.solver.RelTolerance = 1e-5;
fem = fem.solve;