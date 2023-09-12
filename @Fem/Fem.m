classdef Fem < handle
    
    properties (Access = public)
        Mesh;
        materials;
        options;
        solver;
        system;
        triplets;
        topology;

        Dim;           % Dimension (2/3)
        BdBox;          % Bounding box of FEM
        % Mesh;         % Mesh class
        % Material;     % Material class
        % Density;      % Density field (Topology)
        % Node;         % Vertices/Nodes
        % Element;      % Elements
        % Center;       % Element centers
        % NNode;        % Number of nodes
        % NElem;        % Number of elements
        % Topology;     % current save of topology
        % Log;          % Log structure
        % Time;         % Time (solver)
        % Flow;         % auxilary flow function
        % Contraction;  % list of elemental contractions
        % Thickness;    % thickness for 2D fem
        % 
        % Utmp = 0;     % Temporay node diplacements
        % Ftmp = 0;     % Temporay node forces
    end
    
    properties (Access = public, Hidden = true)
        % Node0;         % list of nodes in reference configuration
        % Center0;       % list of element centers in reference configuration
        % Support;       % list of boundary conditions (i.e., supports)
        % Load;          % list of applied forces
        % Tendon;        % list of applied follower forces (i.e., tendons)
        % Spring;        % list of sping attachments to nodes
        % Output;        % list of nodes that are logged
        % Pressure;      % list of edges subjected to pressures
        % PressureCell;  % list of elements undergoing volumetric expansion
        % Contact;       % SDF contact function
        % Source;        
        % 
        % ElemNDof;      % Degrees-of-freedom of elements
        % ElemSort;      % Elements in sorted list
        % ElemRot;       % Rotation mat. each element
        % ElemStr;       % Stretch vector each element
        % Normal;        % Normals of each elements
        % Edge;          % Cell of edges/boundaries
        % Rotation;      % F = RQ: rotational part
        % Stretch;       % F = RQ: stretch part
        % 
        % PrescribedDisplacement = false;
        % VolumetricPressure     = false;
        % TendonForces           = false;
        % GravityRamp            = true;
        % 
        % VonMisesNodal; 
        % sxxNodal; syyNodal; sxyNodal;
        % exxNodal; eyyNodal; exyNodal;
        % fxNodal;  fyNodal;
        % rotNodal;
        % 
        % SPxyNodal;
        % Potential;
        % PotentialG;
        % PotentialF;
        % Kinetic;
        % Gravity;
        % 
        % fContact  = rand(1)*1e-16;
        % fTangent  = rand(1)*1e-16;
        % fInternal = rand(1)*1e-16;
        % fExternal = rand(1)*1e-16;
        % fInput    = rand(1)*1e-16;
        % 
        % Stiffness; 
        % TangentStiffness;
        % SecantStiffness;
        % MassMatrix;
        % DampingMatrix;
        % 
        % Residual = Inf;
        % 
        % TimeEnd      = 1.0;
        % TimeStep     = 0.1;
        % TimeStep0    = 0.1;
        % TimeDelta    = 0;
        % TimeStepMin  = 1e-3;
        % EndIncrement;
        % LoadingFactor;
        % SigmoidFactor;
        % ShowProcess;
        % 
        % ResidualNorm = 1e-3;
        % StressNorm   = 1e-4;
        % DisplaceNorm = 1e-9;
        % Objective    = 1e-4;
        % Constraint   = 1e-4;
        % 
        % Type     = 'PlaneStrain'
        % Solver   = 'mr';
        % SolverId;
        % 
        % IterationMMA = 1;
        % Iteration    = 1;
        % Increment    = 1;
        % Divergence   = 0;
        % Filter       = (1/5); 
        % ShapeFnc;
        % 
        % U        = 0;
        % 
        % dUtmp    = 0;
        % ddUtmp   = 0;
        % Ptmp     = 0;
        % dPtmp    = 0;
        % ddPtmp   = 0;
        % Utmp_    = 0;
        % Ptmp_    = 0;
        % dUtmp_   = 0;
        % ddUtmp_  = 0;
        % Z0       = 0;
        % 
        % PlTmp = 0;
        % 
        % MaxIteration     = 500;
        % MaxIterationMMA  = 50;
        % SolverStart      = false;
        % SolverStartMMA   = false;
        % SolverPlot       = true;
        % SolverPlotType   = 'Svm';
        % QuivPlotContact  = false;
        % Assemble         = false;
        % Convergence      = false; 
        % Nonlinear        = true;
        % BisectLimit      = 15;
        % BisectCounter    = 1;
        % AssembledSystem  = false;
        % PressureLoad     = 0;
        % ParasiticStates  = [];
        % 
        % Linestyle   = '-';
        % Linestyle0  = '-';
        % Colormap    = cmap_turbo;
        % ColormapOpt = cmap_barney(-1);
        % ColorAxis   = [];
        % I3 = eye(3); 
        % O3 = zeros(3);
        % i; j; m; fi; t; e; c; s; p; v; l; k; fb; ed; fb0; ft;
        % 
        % SolverResidual = 1e7;
        % SolverVonMises = 1e7;
        % SolverDisplace = 1e7;
        % 
        % InformationBoolean = false;
        % 
        % Movie        = false;
        % MovieStart   = false;
        % MovieAxis; 
        % MovieCAxis;
        % 
        % VolumeInfill = 0.3;
        % Ersatz       = 1e-3; 
        % Penal        = 1;
        % PenalMax     = 4;
        % PenalStep    = 10;
        % ChangeMax    = 0.03;
        % Beta         = 1.5;
        % FilterRadius = 0.5;
        % SpatialFilter;
        % 
        % OptFactor             = 1;
        % OptimizationSolver    = 'MMA';
        % MaterialInterpolation = 'SIMP';
        % OptimizationProblem   = 'Compliance';
        % xold1; xold2; upp; low; fnorm;
        % zMin = 0; zMax = 1;
        % OutputVector; Change; dFdE;
        % Periodic;
        % 
        % VoidTolerance = 0.05*6; 
        % 
        % Crop;
        % Reflect;
        % NumGridX; NumGridY;
        % Grid; Grideval;
        % Obj = rand(1)*1e-8;
        % Con = rand(1)*1e-8;
        % 
        % TopologyGridX; TopologyGridY; TopologyGridZ;
        % TopologyGapFill;
        % Repeat; ReflectionPlane; CellRepetion;
    end
    
%--------------------------------------------------------------------------
methods  
%---------------------------------------------------------------- Fem Class
function obj = Fem(Mesh,varargin) 

    obj.options   = femoptions;
    obj.topology  = topooptions;
    obj.solver    = solveroptions;
    obj.materials = materialoptions;

    obj.Mesh      = Mesh;
    obj.Dim       = Mesh.get('Dim');
    obj.BdBox     = Mesh.get('BdBox');

    obj.solver.Residual = zeros(obj.Dim*obj.Mesh.NNode,1);
    obj.solver.sol.x    = zeros(obj.Dim*obj.Mesh.NNode,1);
    obj.solver.sol.dx   = zeros(obj.Dim*obj.Mesh.NNode,1);
    obj.solver.sol.ddx  = zeros(obj.Dim*obj.Mesh.NNode,1);
    obj.solver.sol.u    = zeros(obj.Dim*obj.Mesh.NNode,1);

    obj.solver.MaxIteration = 100;

    obj.system.Gravity  = zeros(obj.Dim,1);
    
    if obj.Dim > 2
        obj.options.LineStyle = 'none';
    end
    
    obj.materials      = materiallist;
    obj.topology.sol.x = ones(obj.Mesh.NElem,1);
    
    for ii = 1:2:length(varargin)
        if isprop(obj.options,varargin{ii})
            obj.options.(varargin{ii}) = varargin{ii+1};
        elseif isprop(obj.solver,varargin{ii})
            obj.solver.(varargin{ii}) = varargin{ii+1};
        elseif isprop(obj.topology,varargin{ii})
            obj.topology.(varargin{ii}) = varargin{ii+1};
        else
            obj.(varargin{ii}) = varargin{ii+1};
        end
    end
           
    obj = setupFiniteElement(obj);
end
%---------------------------------------------------------------------- get     
function varargout = get(Fem,varargin)
    if nargin > 1
        varargout{nargin-1,1} = [];
        for ii = 1:length(varargin)
            varargout{ii,1} = Fem.(varargin{ii});
        end
    else
        varargout = Fem.(varargin);
    end
end    
%---------------------------------------------------------------------- set
function Fem = set(Fem,varargin)
Fem = vararginParser(Fem,varargin);
% for ii = 1:2:length(varargin)
%     Fem.(varargin{ii}) = varargin{ii+1};
%     if strcmpi(varargin{ii},'FilterRadius') % regenerate filter
%         Fem.topology.SpatialFilter = GenerateRadialFilter(Fem,varargin{ii+1});
%     end
% end
end
%--------------------------------------------------------------- find nodes
function NodeList = findNodes(Fem,varargin)
    NodeList = findNodeMesh(Fem.Mesh.Node,varargin{1:end});
end
%------------------------------------------------------------ find elements
function ElementList = findElements(Fem,varargin)
    ElementList = findNodeMesh(Fem.Mesh.Center,varargin{1:end});
end
%------------------------------------------------------------ find elements
function NodeList = findEdges(Fem,varargin)
    NodeList = findEdgeMesh(Fem.Mesh,varargin{1:end});
end
%---------------------------------------------------------------- add loads
function Fem = addSpring(Fem,varargin)
 if isa(varargin{1},'char')
    varargin{1} = findNodeMesh(Fem.Mesh.Node,varargin{1});
 end   
 Fem = AddConstraint(Fem,'Spring',varargin{1:end});
end
%--------------------------------------------------------- add displacement
function Fem = addDisplace(Fem,varargin)
 if isa(varargin{1},'char')
     varargin{1} = findNodeMesh(Fem.Mesh.Node,varargin{1});
 end   
 Fem = addDisplaceFem(Fem,varargin{1:end});
end
%-------------------------------------------------------------- add gravity
function Fem = addGravity(Fem,varargin)
 if isempty(varargin)
    if Fem.Dim == 2
        varargin{1} = [0,-9.81e3].';
    else
        varargin{1} = [0,0,-9.81e3].'; 
    end
 end
 Fem = addGravityFem(Fem,[],varargin{1:end});
end
%-------------------------------------------------------- add pressure load
function Fem = addPressure(Fem,varargin)
 if isa(varargin{1},'char')
     varargin{1} = Fem.Mesh.findEdges(varargin{1});
 elseif isempty(varargin{1})
     varargin{1} = Fem.Mesh.findEdges('Hole');
 end   
 Fem = addPressureFem(Fem,varargin{1:end});
end
%------------------------------------------ add myocyte muscle activiation
function Fem = addDilation(Fem,varargin)
    if isa(varargin{1},'char')
        varargin{1} = Fem.Mesh.findElements(varargin{1});
    elseif isempty(varargin{1})
        error('Please assign elements for dilation.');
    end   
    Fem = addDilationFem(Fem,varargin{1:end});
   end
%---------------------------------------------- add volumetric contraction
function Fem = addContraction(Fem,varargin)
    Fem.Contraction = varargin{1};
end
%-------------------------------------- add tendon force or follower force
function Fem = addTendon(Fem,varargin)
    if isa(varargin{1},'char')
        varargin{1} = findNodeMesh(Fem.Mesh.Node,varargin{1});
    end
    Fem = addConstraintFem(Fem,'Tendon',varargin{1:end});
end
%------------------------------------------------ add contact SDF function
function Fem = addContact(Fem,varargin)
 if isa(varargin{1},'Sdf')
     sdf = varargin{1};
     varargin{1} = @(x) sdf.eval(x);
     B1 = box2node(sdf.BdBox); 
     B2 = box2node(Fem.BdBox);
     Fem.BdBox = boxhull([B1;B2],mean(abs(Fem.BdBox))/1e1);
     Fem.system.ContactSDF = sdf;
 end
    
 if numel(varargin) < 2
     varargin{2} = [0,0]; 
 end
 
 Fem = addContactFem(Fem,varargin{1:end});
end
%------------------------------------------------- add output list for LOG
function Fem = addOutput(Fem,varargin)
 if isa(varargin{1},'char')
     varargin{1} = Fem.FindNodes(varargin{1});
 end   
 Fem = AddConstraint(Fem,'Output',varargin{1:end});
end
%----------------------------------------------- add parastic displacements
function Fem = addParastic(Fem,varargin)
[FreeDofs, Ia] = GetFreeDofs(Fem);    
if numel(varargin{1}) == numel(FreeDofs)
    Delta = zeros(Fem.Dim*Fem.NNode,1);
    Delta(Ia) = varargin{1};
    Fem.ParasiticStates = Delta;
else
    Fem.ParasiticStates = varargin{1};
end
end
%---------------------------------------------------- add initial condition
function Fem = addInitialCondition(Fem,varargin)
    if isa(varargin{1},'Fem')
        FEM = varargin{1};
       if ~isempty(FEM.Log)
           Fem.Utmp = FEM.Log.U(end,:).';
       end
    else
       Fem.Utmp = varargin{1}; 
    end
end
%--------------------------------------------------------------- add output
function Fem = addInitialVelocity(Fem,varargin)    
if isempty(Fem.dUtmp)
    Fem.dUtmp = zeros(Fem.NNode*Fem.Dim,1);
end
V0 = varargin{1};
nn = length(Fem.dUtmp)/Fem.Dim;
Fem.dUtmp(1:Fem.Dim:Fem.Dim*nn) = V0(1);
Fem.dUtmp(2:Fem.Dim:Fem.Dim*nn) = V0(2);
if Fem.Dim == 3
    Fem.dUtmp(3:Fem.Dim:Fem.Dim*nn) = V0(3);
end    
end
%---------------------------------- export triangular mesh from iso-surface
function msh = exportMesh(Fem,varargin)
    
    ISO = varargin{1};
    [~,I,UxUy] = showISO(Fem,ISO);
      
    B = UxUy;
    Xscale = (B(2)-B(1))/size(I.CData,2);
    Yscale = (B(4)-B(3))/size(I.CData,1);
    
    simplify_tol = varargin{2};
    
    img = I.CData >= 25;
    img = fliplr(img.');
    
    bnd = bwboundaries(img);
    
    c_cell0 = {};
    c_cell = {};
    
    for ii=1:length(bnd)
        bnd_tmp = bnd{ii};
        assert(all(bnd_tmp(1,:)==bnd_tmp(end,:)),'contour is not closed');
        c_cell0{ii} = bnd_tmp;
    end
    
    for ii=1:length(c_cell0)
        c_tmp = c_cell0{ii};
        c_red = decimatePoly(c_tmp,[simplify_tol, 2],false);
        if (nnz(c_red(:,1))>0)&&(nnz(c_red(:,2))>0)
            c_cell{end+1,1} = [Xscale*c_red(:,1), (Yscale)*c_red(:,2)];
        end
    end
    
    % create the 2d triangulation
    H = varargin{3};
    Tesselation = triangulationCreate(c_cell, H(1), H(2), H(3),'linear');
    
    msh = Mesh(Tesselation.Nodes.',Tesselation.Elements.');
    msh = msh.generate();
end
%----------------------------------------------------------- initial design
function Fem = initialTopology(Fem,varargin)
    if strcmp(varargin{1},'Equidistance')
        Fem.Density = InitialDesign(Fem,{varargin{2},varargin{3}});
    elseif strcmp(varargin{1},'Hole')
        Pc = Fem.Mesh.get('Center');
        d = DistancePointSet(Fem,varargin{2},Pc,varargin{3});
        Z = ones(Fem.NElem,1); Z(d(:,1)) = 0;
        Fem.Density = Z;
    elseif strcmp(varargin{1},'Sdf')
        Pc = Fem.Mesh.get('Center');
        sdf = varargin{2};
        D = sdf(Pc);
        Z = ones(Fem.NElem,1); 
        Z(D(:,end) <= 0) = 0;
        Fem.Density = Z;
    elseif strcmp(varargin{1},'Random')
        Fem.Density = rand(Fem.NElem,1);
    else
        Fem.Density = InitialDesign(Fem,{[1,1],1});
    end
    
    Fem.Density = Fem.SpatialFilter*(1.25*Fem.Density*...
        clamp(Fem.VolumeInfill,0.2,1));
end
%----------------------------------------------------- form smooth topology
function Fem = former(Fem)
Res = 300;
Thickness = 0.2*(Fem.BdBox(2)-Fem.BdBox(1));
Layers = 40;
Patch  = 5;

if nargin < 2, Thickness = 0.2*(Fem.BdBox(2)-Fem.BdBox(1)); end
verts = Fem.Node0; 

fil = full(Fem.SpatialFilter);
V = full(Fem.Mesh.get('NodeToFace'))*fil*Fem.Density;

if Fem.VolumetricPressure
    id = FindElements(Fem,'FloodFill',Fem,Fem.Density);
    Pc = Fem.Mesh.Center;
    S = zeros(Fem.NElem,1);
    d = DistancePointSet(Fem,Pc,Pc(id,:),0.5*Fem.FilterRadius);
    ide = unique(d(:,2));
    
    rhotmp =  fil*Fem.Density;
    S(ide,1) = (rhotmp(ide)>Fem.VoidTolerance)*0.5;
    Fill = zeros(Fem.NElem,1); 
    Fill(id) = 0.5;
    
    E = Fem.Mesh.get('NodeToFace')*...
        fil*(S + Fem.Density*0.01);
    V2 = clamp(Fem.Mesh.get('NodeToFace')*...
        fil*(Fill + Fem.Density-0.1),0,1);
end

if isempty(Fem.Crop)
    B = Fem.BdBox; 
else
    B = Fem.Crop;
end

x = verts(:,1); 
y = verts(:,2);
dx = 0.01*(B(2) - B(1));
dy = 0.01*(B(4) - B(3));
xq = linspace(B(1)+dx,B(2)-dx,Res); 
yq = linspace(B(3)+dy,B(4)-dy,Res);

[xxq,yyq] = meshgrid(xq,yq);
P = griddata(x,y,double(V),xxq,yyq);
vrt = [xxq(:),yyq(:)];
Dist = Fem.Mesh.SDF(vrt); 
Dist = Dist(:,end);

tol = 0.1*sqrt((B(2)-B(1))*(B(4)-B(3))/length(vrt));

P = P(:); 
P(Dist > tol) = 0; 
P = (reshape(P,Res,Res)); 
V = cat(3,xxq,yyq,P);

if Fem.VolumetricPressure
    P = griddata(x,y,double(V2),xxq,yyq);
    Dist = Fem.Mesh.SDF([xxq(:),yyq(:)]); 
    Dist = Dist(:,end);
    
    P = P(:); 
    P(Dist > tol) = 0; 
    P(isnan(P))=0;
    P = (reshape(P,Res,Res));
    
    GapFill = cat(3,xxq,yyq,P);
    
    P = griddata(x,y,double(E),xxq,yyq);
    Dist = Fem.Mesh.SDF([xxq(:),yyq(:)]); 
    Dist = Dist(:,end);
    
    P = P(:);  
    P(isnan(P))=0;
    P = (reshape(P,Res,Res));
    
    EdgeFill = cat(3,xxq,yyq,P);
end

SDF0 = V(:,:,3);
SDF0 = GaussianFilter(SDF0,5);
SDF  = repmat(SDF0,[1 1 Layers]);

if Fem.VolumetricPressure
    V       = GapFill;
    SDF2    = V(:,:,3);
    SDF2    = GaussianFilter(SDF2,10);
    ZFiller = repmat(SDF2,[1 1 Patch]);
    
    for ii = 1:Patch
        ZFiller(:,:,ii) = lerp(SDF2+0.01*ii,SDF0,...
                               cos((1/Patch)*pi));
    end
end

% zero padding
SDF(:,:,1) = 0; SDF(:,:,end) = 0;
SDF(1,:,:) = 0; SDF(end,:,:) = 0;
SDF(:,1,:) = 0; SDF(:,end,:) = 0;

Fem.Topology = SDF;

zq = linspace(0,Thickness,Layers);
[X,Y,Z] = meshgrid(xq,yq,zq);

Fem.TopologyGridX = X; 
Fem.TopologyGridY = Y; 
Fem.TopologyGridZ = Z; 

end
%--------------------------------------------------------- show iso-surface
function [Fem,I,UxUy] = showISO(Fem,varargin)
V = Fem.Topology;
X = Fem.TopologyGridX; 
Y = Fem.TopologyGridY; 
a = size(V);

if nargin < 3
    depth = floor(a(3)/2);
else
    depth = a(3)*varargin{2}; 
end

X    = X(:,:,depth);
Y    = Y(:,:,depth);
SDF0 = V(:,:,depth);

SDF = GenerateCell(Fem,SDF0);

if ~isempty(Fem.CellRepetion)
    Rpt = Fem.CellRepetion;
    SDF = repmat(SDF,Rpt(2),Rpt(1));
else
    Rpt = [1, 1];
end

if ~isempty(Fem.ReflectionPlane)
    Rp = Fem.ReflectionPlane;
    
    if Rp(1) == 1
        Xf = flip(X) + (max(max(max(X))) - min(min(min(X))));
        X  = vertcat(Xf,X);
    end
    
    if Rp(1) == -1
        Xf = flip(X) - (max(max(max(X))) - min(min(min(X))));
        X  = horzcat(X,Xf);
    end
    
    if Rp(2) == -1
        Yf = flip(Y) - (max(max(max(Y))) - min(min(min(Y))));
        Y  = horzcat(Y,Yf);
    end
    
    if Rp(2) == 1
        Yf = flip(Y) + (max(max(max(Y))) - min(min(min(Y))));
        Y  = horzcat(Yf,Y);
    end
end

scaleX = 1; scaleY = 1;
if ~isempty(Fem.Repeat)
    
    SDF0     = SDF;
    Instruct = Fem.Repeat;
    
    for ii = 1:length(Instruct)
        
        if Instruct(ii) == 1
            scaleX = scaleX + 1;
            SDF    = cat(2,SDF,SDF0);
        end
        
        if Instruct(ii) == 2
            scaleY = scaleY*2;
            SDF    = cat(1,SDF,SDF);
        end
    end
end

cla;
Uxx  = scaleX*Rpt(1)*[min(min(min(X))) max(max(max(X)))];
Uyy  = scaleY*Rpt(2)*[min(min(min(Y))) max(max(max(Y)))];
UxUy = [0, Uxx(2)-Uxx(1), 0, Uyy(2)-Uyy(1)];

I = GaussianFilter((SDF >= varargin{1})*255,3);

I = image(rescale(Uxx),...
    ((max(Uyy)-min(Uyy))/(max(Uxx) - min(Uxx)))*rescale(Uyy),I);

axis equal; axis off; 
colormap(Fem.ColormapOpt); 
caxis([0 1]);
background('w');
end


end
end