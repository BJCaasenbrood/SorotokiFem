clr;

sdf = sCircle(60) - sCircle(57);
msh = Mesh(sdf,'NElem',600);
msh = msh.generate();

fem = Fem(msh,'TimeStep',1/500);
fem = fem.addMaterial(NeoHookean(1.0,0.33));
fem = fem.addPressure('allhole',100 * 1e-3);

fem = fem.addSupport('Bottom',[1,1]);

fem = solveQuasiStaticFem(fem);