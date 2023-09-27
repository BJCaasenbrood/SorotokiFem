function [Fnc, Ftc, Knc, Ktc] = assembleContactForces(Fem)

omegaN = Fem.materials.Material{1}.getContactReaction;
omegaT = Fem.materials.Material{1}.getContactTangentReaction;
mu     = Fem.materials.Material{1}.getContactFriction;

Knc = sparse(Fem.Mesh.NNode*Fem.Dim,Fem.Mesh.NNode*Fem.Dim);
Ktc = sparse(Fem.Mesh.NNode*Fem.Dim,Fem.Mesh.NNode*Fem.Dim);
Fnc = sparse(Fem.Mesh.NNode*Fem.Dim,1);  
Ftc = sparse(Fem.Mesh.NNode*Fem.Dim,1);  
SDF = Fem.system.Contact{1};

Ia = Fem.system.ContactMesh.NodeId;
qa = Fem.system.ContactMesh.qa;

Y0   = Fem.Mesh.Node(Ia,:);
[U]  = meshfield(Fem,Fem.solver.sol.x);
[dU] = meshfield(Fem,Fem.solver.sol.dx);

Y   = Y0 + U(Ia,:);
eps = Fem.solver.TimeStep * 10;

d = SDF(Y);
Intersect = find((d(:,end))<eps);

if ~isempty(Intersect)
    I = Intersect;
    
    n1 = (SDF(Y(I,:)+repmat([eps,0],size(Y(I,:),1),1))-d(I,end))/eps;
    n2 = (SDF(Y(I,:)+repmat([0,eps],size(Y(I,:),1),1))-d(I,end))/eps;

    % normal vector, tangent vector
    N(:,1) =  n1(:,end);
    N(:,2) =  n2(:,end);
    T(:,2) =  n1(:,end);
    T(:,1) = -n2(:,end);

    V0 = dU(Ia(I),:);
    % Vn = V0./(vecnorm(V0.').');

    % friction penalty function (i.e, expected distance traveled)
    gt = (dot(V0.',T.').') * Fem.solver.TimeStep;

    % contact penalty function
    dist = abs(d(I,end));
    % dhat = clamp(abs(dot(V0.', N.').') * Fem.solver.TimeStep,0,1);
    % ipcType =  dist < dhat;

    gn = -dist;
    % gn = ipcType .* ( - (dist - dhat).^2 .* log(dist./dhat) ) ./ (dist/Fem.solver.TimeStep);

    Ux = gn.*n1(:,end);
    Uy = gn.*n2(:,end);

    % stick-slip boolean vector
    cType = ~(abs(omegaT * gt) >= abs(mu * omegaN * gn));

    Fcont = -omegaN * [Ux, Uy];
    FfricStick = -omegaT*gt.*T;
    FfricSlip = mu*omegaN*sign(gt).*gn.*T;

    % RF = -omegaN * [Ux, Uy];
    % RF = (vecnorm(RF.').');

    Ix = 2*Ia(I(1:size(Y(I,:),1),1))-1;
    Iy = 2*Ia(I(1:size(Y(I,:),1),1));

    Fnc(Ix,1) = Fcont(:,1);
    Fnc(Iy,1) = Fcont(:,2);
    Ftc(Ix,1) = FfricStick(:,1) .* cType +  FfricSlip(:,1) .* (1-cType);
    Ftc(Iy,1) = FfricStick(:,2) .* cType +  FfricSlip(:,2) .* (1-cType);
    
    index = 0;
    NDof  = Fem.Dim;
    I = zeros(numel(Ix)*NDof^2,1);
    J = zeros(numel(Ix)*NDof^2,1);
    Ke = zeros(numel(Ix)*NDof^2,1);
    Te = zeros(numel(Ix)*NDof^2,1);

    for ii = 1:numel(Ix)
        
        KK = omegaN*N(ii,:) .' * N(ii,:);
        TT = cType(ii)*(omegaT*T(ii,:).'*T(ii,:)) + ...
            (1-cType(ii))*(mu*omegaN*sign(gt(ii))*T(ii,:).'*N(ii,:));

        eDof = [Ix(ii); Iy(ii)];
        II = repmat(eDof,1,Fem.Dim);
        JJ = II';

        ind = index+1:index+NDof^2;
        I(ind) = II(:);
        J(ind) = JJ(:);
        Ke(ind) = KK(:);
        Te(ind) = TT(:);

        index = index + NDof^2;
    end

    Knc = sparse(I,J,Ke, ...
        Fem.Mesh.NNode*Fem.Dim,Fem.Mesh.NNode*Fem.Dim);

    Ktc = omegaT*sparse(I,J,Te, ...
        Fem.Mesh.NNode*Fem.Dim,Fem.Mesh.NNode*Fem.Dim);
end

end
