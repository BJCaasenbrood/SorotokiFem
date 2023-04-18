clr;
W = 150;
H = 12;

%% generate mesh
msh = Mesh('Crawler.png','BdBox',[0,W,0,H],'NElem',750);
msh = msh.generate();
msh.show();

% %% FEM
fem = Fem(msh,'TimeStep',1/1250,'TimeHorizon',6);

fem.addMaterial(NeoHookean(0.2,0.4));
fem.addMaterial(NeoHookean(1,0.4));

inExtLayer = fem.findElements('Box',[0,150,0,2]);
fem = fem.setMaterial(inExtLayer,2);

% fem.Material.Zeta = .75;
% fem.Material.Rr   = .3;
% fem.Material.Cfr  = .5;

fem = fem.addGravity();
fem = fem.addContact(SDF(W));

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
    %warning on;
    %axis([-50 120 -100 20]);
end

function y = pset(t,id)
w  = 2*pi;
P0 = 15 * 1e-3;
t  = mod(t,.75);
y  = P0*tsin(w*t - (id-1)*pi/5).*sigmoid(3*w*t);
end
    

% %% make inextensible bottom layer
% id = fem.FindElements('Box',[0,150,0,2]);
% fem.Density(id) = 5;

% %%
% fem = fem.addPressure(fem.FindEdges('BoxHole',[0 50 0 15]),...
%    @(x) pset(x.Time,1));
% fem = fem.addPressure(fem.FindEdges('BoxHole',[50 100 0 15]),...
%    @(x) 0.75*pset(x.Time,2));
% fem = fem.addPressure(fem.FindEdges('BoxHole',[100 150 0 15]),...
%    @(x) pset(x.Time,3));

% %% simulate
% fem = fem.simulate();

% function y = pset(t,id)
% w  = 2*pi;
% P0 = 15 * kpa;
% t  = mod(t,.75);
% y  = P0*tsin(w*t - (id-1)*pi/5).*sigmoid(3*w*t);
% end

function y = SDF(W)
    y = sLine(1.5*W,-0.1*W,.01,.01);
end