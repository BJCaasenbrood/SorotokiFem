function [E, dEdy, V, dVdy] = materialFieldFem(Fem, ForceFilterOff)

    if nargin < 2
        ForceFilterOff = false;
    end

    if ~ForceFilterOff 
      y = Fem.topology.SpatialFilter * ...
            Fem.topology.sol.x;
    else
        y = Fem.topology.sol.x;
    end

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
