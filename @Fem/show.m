function h = show(Fem,varargin)

Nds         = Fem.Mesh.Node + meshfield(Fem, Fem.solver.sol.x);
FaceMatrix  = Fem.Mesh.geometry.ElemMat;
BoundMatrix = Fem.Mesh.geometry.Boundary;
warning off;

h{1} = patch('Faces', FaceMatrix, 'Vertices', Nds, ...
    'LineStyle', 'none', 'Linewidth', 1.5, 'FaceColor',...
    [0.95,0.95,0.957]);

h{3} = patch('Faces',BoundMatrix,'Vertices',Nds,...
    'LineStyle','-','Linewidth',1.5,'EdgeColor','k');

h{2} = patch('Faces', FaceMatrix, 'Vertices', Nds, ...
    'Facecolor', Fem.options.Color, 'LineStyle',             ...
    Fem.options.LineStyle,'Linewidth', 1.5, 'FaceAlpha', 1,...
    'EdgeColor', 'k');

axis equal;
axis off; 
hold on;
warning on;
end
