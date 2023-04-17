%--------------------------------------------------------------------- show
function h = showMaterialsFem(Fem)

Z = Fem.materials.MatElem + 1;


h{3} = [];

FaceMatrix  = Fem.Mesh.geometry.ElemMat;
BoundMatrix = Fem.Mesh.geometry.Boundary;
warning off;
if Fem.Dim < 3
    h{1} = patch('Faces', FaceMatrix, 'Vertices', Fem.Mesh.Node, ...
        'LineStyle', 'none', 'Linewidth', 1.5, 'FaceColor', 'flat');
end

h{3} = patch('Faces',BoundMatrix,'Vertices',Fem.Mesh.Node,...
    'LineStyle','-','Linewidth',1.5,'EdgeColor','k');

h{2} = patch('Faces', FaceMatrix, 'Vertices', Fem.Mesh.Node, ...
             'FaceVertexCData', Z, 'Facecolor', 'flat', 'LineStyle',             ...
             Fem.options.LineStyle,'Linewidth', 1.5, 'FaceAlpha', 1,...
             'EdgeColor', 'k');

axis equal;     
axis off; hold on;
colormap(Fem.materials.ColorMap); clim([0,7]);
warning on;

end
