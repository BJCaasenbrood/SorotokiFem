function h = showMaterialsFem(Fem)

h = Fem.Mesh;
Z = Fem.materials.MatElem;

h = {1};

Nds         = Fem.Mesh.Node + meshfield(Fem, Fem.solver.sol.x);
FaceMatrix  = Fem.Mesh.geometry.ElemMat;
BoundMatrix = Fem.Mesh.geometry.Boundary;
warning off;

if Fem.Dim < 3
    h{1} = patch('Faces', FaceMatrix, 'Vertices', Nds, ...
        'LineStyle', 'none', 'Linewidth', 1.5, 'FaceColor',...
        [0.95,0.95,0.957]);
end

h{3} = patch('Faces',BoundMatrix,'Vertices',Nds,...
    'LineStyle','-','Linewidth',1.5,'EdgeColor','k');

h{2} = patch('Faces', FaceMatrix, 'Vertices', Nds, ...
             'FaceVertexCData', Z, 'Facecolor', 'flat', 'LineStyle',             ...
             Fem.options.LineStyle,'Linewidth', 1.5, 'FaceAlpha', 1,...
             'EdgeColor', 'k');

axis equal;     
axis off; hold on;
colormap(Fem.materials.ColorMap); 
clim([0.5,6.5]);
colorbar;
warning on;

end
