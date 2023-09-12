% 10-Apr-2023 11:50:28
% Auto-generated test script

% Initialize the test suite
% Add test cases here
clr;
s = sRectangle(-5,5,0,40);
msh = Mesh(s,'Quads',[10,30]);
msh = msh.generate();

E = msh.findElements('Box',[-3,3,4,36]);
msh = msh.removeElements(E);

fem = Fem(msh,'TimeStep',1/10);
fem = fem.addMaterial(NeoHookean);
fem = fem.addMaterial(NeoHookean(100,0.3));
fem = fem.addMaterial(NeoHookean(10,0.3));

fem = fem.setMaterial('left',2);
fem = fem.setMaterial('right',2);
fem = fem.setMaterial(msh.findElements('box',[-5,5,0,4]),3);
fem = fem.setMaterial(msh.findElements('box',[-5,5,36,40]),3);

I = fem.findEdges('hole',2);
fem = fem.addPressure(I,15 * 1e-3);
fem = fem.addSupport('top',[1,1]);
qa = fem.system.Ia;

fem = fem.set('isLog',false);

% fem = solveQuasiStaticFem(fem);
% showVonMisesFem(fem);