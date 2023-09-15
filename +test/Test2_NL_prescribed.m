% 10-Apr-2023 08:10:35
% Auto-generated test script

% Initialize the test suite
% Add test cases here

sdf = sRectangle(0, 100, 0, 50);
msh = Mesh(sdf,'Quads',[30,3]);
msh = msh.generate();

fem = Fem(msh,'TimeStep',1/10);
fem = fem.addMaterial(Yeoh([1]));
fem = fem.addSupport('left', [1, 1]);
fem = fem.addDisplace('right', [250,0]);

fem.options.isNonlinear = true;
fem.solver.MaxIteration = 20;

fem = fem.solve;
plt(fem);

function plt(Fem)
    showVonMisesFem(Fem);
    axis tight;
end
