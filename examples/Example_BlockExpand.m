sdf = sRectangle(0,1,0,1);
msh = Mesh(sdf,'Quads',[9,9]);
msh = msh.generate();

msh.show();

fem = Fem(msh,'TimeStep',1/30);
fem = fem.addSupport('bottom',[1,1]);
fem = fem.addMaterial(Yeoh(0.1));

id = fem.findElements('Box',[0.333 0.666 0.333 0.666]);
fem = fem.addDilation(id,-0.5);

fem = fem.solve;