function h = showVonMisesFem(Fem,varargin)

Z = full(sparse(Fem.triplets.l,1,Fem.triplets.vj));

Nds         = Fem.Mesh.Node + meshfield(Fem, Fem.solver.sol.x);
FaceMatrix  = Fem.Mesh.geometry.ElemMat;
BoundMatrix = Fem.Mesh.geometry.Boundary;
warning off;

if isempty(varargin)
    h = {};

    % if Fem.Dim < 3
    h{1} = patch('Faces', FaceMatrix, 'Vertices', Nds, ...
        'LineStyle', 'none', 'Linewidth', 1.5, 'FaceColor',...
        [0.95,0.95,0.957]);
    % end
    %
    h{3} = patch('Faces',BoundMatrix,'Vertices',Nds,...
        'LineStyle','-','Linewidth',1.5,'EdgeColor','k');

    h{2} = patch('Faces', FaceMatrix, 'Vertices', Nds, ...
        'FaceVertexCData', Z, 'Facecolor', 'interp', 'LineStyle',             ...
        Fem.options.LineStyle,'Linewidth', 1.5, 'FaceAlpha', 1,...
        'EdgeColor', 'k');

    axis equal;
    axis off; hold on;
    colormap(flipud(Fem.options.ColorMap));
    caxis([-.5,.5]);
    warning on;
else    

    set(h{2},'Vertices',Nds);
    set(h{2},'FaceVertexCData',Z);

end
end

