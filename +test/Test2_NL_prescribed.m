% 10-Apr-2023 08:10:35
% Auto-generated test script

% Initialize the test suite
% Add test cases here

clr;
sdf = sRectangle(0, 100, 0, 50);
msh = Mesh(sdf,'Quads',[30,3],'MaxIteration',50);
msh = msh.generate();

fem = Fem(msh,'TimeStep',1/10);
fem = fem.addMaterial(NeoHookean(1,0.4));
fem = fem.addSupport('left', [1, 1]);
fem = fem.addDisplace('right', [250,0]);

fem.options.isNonlinear = true;

fem = solveQuasiStaticFem(fem);
plt(fem);

function plt(Fem)
    showVonMisesFem(Fem);
    axis tight;
end
