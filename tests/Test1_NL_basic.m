% 10-Apr-2023 08:10:35
% Auto-generated test script

% Initialize the test suite
% Add test cases here
clr;

sdf = sRectangle(0, 100, 0, 5);
msh = Mesh(sdf,'NElem',120,'MaxIteration',150);
msh = msh.generate();

fem = Fem(msh,'TimeStep',0.1);
fem = fem.addMaterial(Ecoflex0030);
fem = fem.addSupport('left', [1, 1]);
fem = fem.addSupport('right',[1, 1]);
fem = fem.addLoad('bottom', [0, 5e-4]);

fem = solveQuasiStaticFem(fem);
