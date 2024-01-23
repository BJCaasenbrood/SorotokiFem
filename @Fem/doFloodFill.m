function out = doFloodFill(Fem)

    Nds = Fem.Mesh.Node;
    rho = Fem.topology.sol.x;
    adj = Fem.Mesh.geometry.Adjecency;

    eps = 0.1;
    
    if isfield(Fem.options,'VoidTolerance')
        eps = Fem.options.VoidTolerance;
    end
    
    blclist = find(rho > eps);
    adj(blclist,:) = 0;
    adj(:,blclist) = 0;
    
    % injectionPoint = Fem.get('PressureCell');
    injectionPoint = Fem.system.Myocyte{1};
    
    x   = FloodFillWhileLoop(adj, injectionPoint(1));
    id  = 1:length(Nds);
    out = id( logical(x) );

end

function list = FloodFillWhileLoop(adj,v0)
[n,~] = size(adj);
list = zeros(n,1);
list(v0) = 1;
neigh = (adj(:,v0) & ~list);
stack = find(neigh);

    while(~isempty(stack))
        v0 = stack(end);
        stack(end) = [];
        list(v0)=1;
        neigh=(adj(:,v0) & ~list); 
        stack = unique( [stack; find(neigh)]);  
    end

end
