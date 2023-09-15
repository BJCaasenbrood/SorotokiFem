sdf = sRectangle(0,1,0,1);
msh = Mesh(sdf,'Quads',[9,9]);
msh = msh.generate();

msh.show();

fem = Fem(msh,'TimeStep',1/30);
fem = fem.addSupport('Left',[1,1]);
% fem = fem.addMaterial(NeoHookean(0.1,0.499));
fem = fem.addMaterial(Yeoh(0.1));

fem = fem.addDisplace('right',[5,0]);

fem = fem.solve;