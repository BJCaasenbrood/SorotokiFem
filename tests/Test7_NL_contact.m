% 17-Apr-2023 17:09:43
% Auto-generated test script

% Initialize the test suite
% Add test cases here

% mesh
con = SDF(15);
sdf = sCircle(5,[18,25]) - sCircle(3,[18,25]);
msh = Mesh(sdf,'NElem',55,'MaxIteration',500);
msh = msh.generate();

msh.BdBox = [-20,30,-10,30];

fem = Fem(msh,'TimeStep',1/850,'TimeHorizon',.1);
fem = fem.addMaterial(NeoHookean(0.01,0.4));

fem = fem.addGravity();
fem = fem.addContact(con);

fem.options.Display = @plt;
fem = solveDynamicFem(fem);

function plt(Fem)
    cla;
    showVonMisesFem(Fem);
    showContactFem(Fem);
    xlim([-20,30]);
    ylim([-10,30]);
end

function D = SDF(W)
    S1 = sLine(4,-W,-10,-5);
    S2 = sLine(W,4,0,-10);
    S3 = sLine(-W,-W,-1,1);
    D = S1 + S2 + S3;
end