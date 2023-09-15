% 11-Apr-2023 11:40:37
% Auto-generated test script

% Initialize the test suite
% Add test cases here

t = 1.5;
h = 7.5;
l = 75;
w = h+3*t;

s = sRectangle(0,3*t + h,0,l);
msh = Mesh(s,'Quads',[6,70]);

msh = msh.generate();

NBellow = 8;
X0 = 10;
W  = 5;
dX = (l - 2*X0 - NBellow*W)/(NBellow-1);

for ii = 1:8
    Pc = msh.Center;
    sR = sRectangle(t,t+h,X0 + (W+dX)*(ii-1), X0 + (W+dX)*(ii-1) + W);
    [~,E] = sR.intersect(Pc);
    msh = msh.removeElements(E);
end

fem = Fem(msh,'TimeStep',1/60);
fem = fem.addMaterial(preset.material.Ecoflex0030(20));
fem = fem.addMaterial(NeoHookean(1,0.3));

fem = fem.setMaterial(msh.findElements('box',[w-2*t,w,0,l]),2);

showMaterialsFem(fem); axis on;

fem = fem.addPressure(fem.findEdges('allhole'), 25 * 1e-3);
fem = fem.addSupport('bottom',[1,1]);

fem = fem.solve;
