clr;
W = 150;
H = 12;

%% generate mesh
msh = Mesh('Crawler.png','BdBox',[0,W,0,H],'NElem',500);
msh = msh.generate();
msh.show();

% %% FEM
fem = Fem(msh,'TimeStep',1/350,'TimeHorizon',6);

mat = NeoHookean(0.2,0.4);
mat.contact.TangentFriction = 0.5;

fem.addMaterial(mat);
fem.addMaterial(NeoHookean(1,0.4));

BotLayer = fem.findElements('Box',[0,150,0,2]);
fem = fem.setMaterial(BotLayer,2);

fem = fem.addGravity();
fem = fem.addContact(SDF);

fem.options.LineStyle = 'none';
fem.options.Display   = @plt;

fem = fem.addPressure(fem.findEdges('BoxHole',[0 50 0 15]),...
   @(t) pset(t,1));
fem = fem.addPressure(fem.findEdges('BoxHole',[50 100 0 15]),...
   @(t) 0.75*pset(t,2));
fem = fem.addPressure(fem.findEdges('BoxHole',[100 150 0 15]),...
   @(t) pset(t,3));

fem = solveDynamicFem(fem);

function plt(Fem)
    clf;
    %warning off;
    showVonMisesFem(Fem);
    showContactFem(Fem);
    axis([-20 260 -5 50]);
    axis on;
end

function y = pset(t,id)
w  = 2*pi;
P0 = 13 * 1e-3;
t  = mod(t,.75);
y  = P0*tsin(w*t - (id-1)*pi/5).*sigmoid(3*w*t);
end
    

function y = SDF
   y = sLine(250,-10,.01,.01);
end