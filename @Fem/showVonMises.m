function h = showVonMises(Fem,varargin)

% compose the stress values
Z = full(sparse(Fem.triplets.l,1,Fem.triplets.s(:,1)));

Nds         = Fem.Mesh.Node + meshfield(Fem, Fem.solver.sol.x);
FaceMatrix  = Fem.Mesh.geometry.ElemMat;
BoundMatrix = Fem.Mesh.geometry.Boundary;

if length(Z) == Fem.Mesh.NNode
    Z = Smoothing(Fem,Z,3); 
end

warning off;

h{1} = patch('Faces', FaceMatrix, 'Vertices', Nds, ...
    'LineStyle', 'none', 'Linewidth', 1.5, 'FaceColor',...
    [0.95,0.95,0.957]);

h{3} = patch('Faces',BoundMatrix,'Vertices',Nds,...
    'LineStyle','-','Linewidth',1.5,'EdgeColor','k');

h{2} = patch('Faces', FaceMatrix, 'Vertices', Nds, ...
    'FaceVertexCData', Z, 'Facecolor', 'interp', 'LineStyle',             ...
    Fem.options.LineStyle,'Linewidth', 1.5, 'FaceAlpha', 1,...
    'EdgeColor', 'k');

axis off; hold on;
colormap(Fem.options.ColorMap);
warning on;
end

%% mesh smoothing
function f = Smoothing(Fem,f,naver)
    W = NodeAdjecency(Fem.Mesh.Element);
    n = length(W); 
    D = spdiags(full(sum(W,2).^(-1)),0,n,n);
    W = D*W;
    
    for ii=1:naver
        f = W*double(f(:));
    end
end

%% GENERATE ADJECENCY MATRIX FROM FACES
function A = NodeAdjecency(face)
n = max(cellfun(@max,face));      
A = sparse(n,n);
null = sparse(n,n);

for ii = 1:length(face)
    poly = uint16(face{ii});
    [i,j] = PermutationSet(numel(poly));
    B = null;
    B(poly(i),poly(j)) = 1;
    A = A + B;
end

A = double(A>0);
end

%% PERMUTATION SET
function [i,j] = PermutationSet(n)
i = transpose(kron(1:n,ones(1,n)));
j = transpose(kron(ones(1,n),1:n));
end
