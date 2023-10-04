function h = showInfillFem(Fem,varargin)

[~,~,Z]   = materialField(Fem);

Nds         = Fem.Mesh.Node;
FaceMatrix  = Fem.Mesh.geometry.ElemMat;
BoundMatrix = Fem.Mesh.geometry.Boundary;
warning off;

h = {};

h{1} = patch('Faces', FaceMatrix, 'Vertices', Nds, ...
    'LineStyle', 'none', 'Linewidth', 1.5, 'FaceColor',...
    [0.95,0.95,0.957]);

h{3} = patch('Faces',BoundMatrix,'Vertices',Nds,...
    'LineStyle','-','Linewidth',1.5,'EdgeColor','k');

h{2} = patch('Faces', FaceMatrix, 'Vertices', Nds, ...
    'FaceVertexCData', Z, 'Facecolor', 'flat', 'LineStyle',             ...
    Fem.options.LineStyle,'Linewidth', 1.5, 'FaceAlpha', 1,...
    'EdgeColor', 'k');

axis equal;
axis off; hold on;
colormap(Fem.topology.ColorMap);
clim([Fem.topology.Ersatz,1]);
warning on;

end

function [E, dEdy, V, dVdy] = materialField(Fem)

    y = Fem.topology.SpatialFilter * Fem.topology.sol.x;
    eps = Fem.topology.Ersatz;

    switch (Fem.topology.Interpolation)
        case ('SIMP')
            penal = clamp(Fem.topology.Penal, 1, Fem.topology.MaxPenal);
            E = eps + (1 - eps) * y .^ penal;
            V = y;
            dEdy = (1 - eps) * penal * y .^ (penal - 1);
            dVdy = ones(size(y, 1), 1);
        case ('SIMP-H')
            penal = clamp(Fem.topology.Penal, 1, Fem.topology.MaxPenal);
            beta = Fem.topology.Beta;
            h = 1 - exp(-beta * y) + y * exp(-beta);
            E = eps + (1 - eps) * h .^ penal;
            V = h;
            dhdy = beta * exp(-beta * y) + exp(-beta);
            dEdy = (1 - eps) * penal * h .^ (penal - 1) .* dhdy;
            dVdy = dhdy;
    end
end


