classdef Fem < handle
    
    properties (Access = public)
        Mesh;
        log;
        materials;
        options;
        solver;
        system;
        triplets;
        topology;

        Dim;   % Dimension (2/3)
        BdBox; % Bounding box of FEM
    end
    
    properties (Access = public, Hidden = true)

    end
    
%--------------------------------------------------------------------------
methods  
%---------------------------------------------------------------- Fem Class
function obj = Fem(Mesh,varargin) 

    obj.log       = Log;
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
%--------------------------------------------------------- add displacement
function Fem = addDisplace(Fem,varargin)
 if isa(varargin{1},'char')
     varargin{1} = findNodeMesh(Fem.Mesh.Node,varargin{1});
 end   
 Fem = addDisplaceFem(Fem,varargin{1:end});
end
%---------------------------------------------- add volumetric contraction
% function Fem = addContraction(Fem,varargin)
%     Fem.Contraction = varargin{1};
% end
%-------------------------------------- add tendon force or follower force
function Fem = addTendon(Fem,varargin)
    if isa(varargin{1},'char')
        varargin{1} = findNodeMesh(Fem.Mesh.Node,varargin{1});
    end
    Fem = addConstraintFem(Fem,'Tendon',varargin{1:end});
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

end
end