clr;
warning off;

sdf = sRectangle(0, 100, 0, 3);
msh = Mesh(sdf,'Quads',[50,1],'MaxIteration',150);
msh = msh.generate();

fem = Fem(msh,'TimeStep',1/1250);
fem = fem.addMaterial(NeoHookean(.3,0.49));
fem = fem.addSupport('left', [1, 1]);
fem = fem.addContact(sCircle(7,[40,-25]));
fem = fem.addGravity();

fem.options.LineStyle   = 'none';
%fem.options.isNonlinear = false;
fem.options.Display = @plt;

fem = solveDynamicFem(fem);

function plt(Fem)
    clf;
    %warning off;
    showVonMisesFem(Fem);
    showContactFem(Fem);
    %warning on;
    %axis([-50 120 -100 20]);
end