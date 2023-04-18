function F = addBasicLoadsFem(Fem)

    F = sparse(Fem.Mesh.NNode*Fem.Dim,1);   
    Nds = Fem.Mesh.Node;
    X   = Fem.solver.sol.x;

    Nds(:, 1) = Nds(:, 1) + X(1:2:end - 1, 1);
    Nds(:, 2) = Nds(:, 2) + X(2:2:end, 1);

    for kk = 1:size(Fem.system.Pressure,1)

        EdgeList = Fem.system.Pressure{kk, 1};
        Pload = Fem.system.Pressure{kk, 2}(Fem.solver.Time);

        % if ~isa(Fem.system.Pressure{kk, 2}, 'function_handle')
        %     Pload = mean(Fem.system.Pressure{kk, 2});
        % else
        %     Pload = Fem.system.Pressure{kk, 2}(Fem);
        % end

        for ii = 1:length(EdgeList)
            NodeID = EdgeList{ii};
            V = Nds(NodeID, :);

            % [~,N] = gradient(V);
            %     N = N./sqrt((sum((N.^2),2)) + 1e-3 );
            %     Ny = N(:,1);
            %     Nx = N(:,2);

            dsx = diff(V(:, 1));
            dsy = diff(V(:, 2));
            dl = sqrt(dsx .^ 2 + dsy .^ 2);
            Nx = -dsy ./ dl;
            Ny = dsx ./ dl;

            S = [NodeID(1:end - 1), NodeID(2:end)].';

            for jj = 1:length(NodeID) - 1
                F(Fem.Dim * S(1, jj) - 1, 1) = F(Fem.Dim * S(1, jj) - 1, 1) ...
                    + 0.5 * Pload * Nx(jj) * dl(jj);
                F(Fem.Dim * S(2, jj) - 1, 1) = F(Fem.Dim * S(2, jj) - 1, 1) ...
                    + 0.5 * Pload * Nx(jj) * dl(jj);
                F(Fem.Dim * S(1, jj), 1) = F(Fem.Dim * S(1, jj), 1) ...
                    + 0.5 * Pload * Ny(jj) * dl(jj);
                F(Fem.Dim * S(2, jj), 1) = F(Fem.Dim * S(2, jj), 1) ...
                    + 0.5 * Pload * Ny(jj) * dl(jj);

            end

        end

    end

end
