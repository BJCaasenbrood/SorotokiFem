function [g, dgdE, dgdV] = ConstraintFunction(Fem)

    Vd = Fem.topology.VolumeInfill;
    V  = Fem.topology.SpatialFilter * Fem.topology.sol.x;
    A  = Fem.Mesh.get('Area');

    g = sum(A .* V) / (sum(A) * Vd) - 1;
    dgdE = zeros(size(A));
    dgdV = A / (sum(A) * Vd);

end
