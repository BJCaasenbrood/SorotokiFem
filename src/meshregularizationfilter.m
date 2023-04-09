function P = meshregularizationfilter(Mesh, Radius, varargin)

Periodic = [];
if strcmpi(varargin{1},'periodic')
    Periodic = varargin{2};
end

if Radius < 0
    P = speye(Mesh.NElem); 
    return;
end

PS1 = Mesh.get('Center');
PS2 = Mesh.get('Center');

d = distancepointset(PS1,PS2,Radius,                                         ...
        'boundingBox',Mesh.BdBox,                                            ...
        'periodicBoundary',Periodic);
    
if ~isempty(d)
    P = sparse(d(:,1),d(:,2),1-d(:,3)/Radius);
    P = spdiags(1./sum(P,2),0,Mesh.NElem,Mesh.NElem)*P;
else
    P = 1;
end
end

