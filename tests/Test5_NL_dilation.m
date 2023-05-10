% 10-Apr-2023 11:50:28
% Auto-generated test script

% Initialize the test suite
% Add test cases here

clf;
s = sRectangle(-5,5,0,40);
msh = Mesh(s,'Quads',[5,10]);
msh = msh.generate();

E = msh.findElements('Box',[-3,3,4,36]);

fem = Fem(msh,'TimeStep',1/130);
fem = fem.addMaterial(NeoHookean);
fem = fem.addMaterial(NeoHookean(0.01,0.49));
fem = fem.setMaterial(E,2);
fem = fem.addSupport('top',[1,1]);
fem = fem.addSupport('bottom',[0,1]);
fem = fem.addDilation(E,-.4);

fem = solveQuasiStaticFem(fem);
showVonMisesFem(fem);
