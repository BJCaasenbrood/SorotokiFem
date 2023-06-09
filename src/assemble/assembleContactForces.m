function [Fnc, Ftc, Knc, Ktc] = assembleContactForces(Fem)

omegaN = Fem.materials.Material{1}.getContactReaction;
omegaT = Fem.materials.Material{1}.getContactFriction;

Knc = sparse(Fem.Mesh.NNode*Fem.Dim,Fem.Mesh.NNode*Fem.Dim);
Ktc = sparse(Fem.Mesh.NNode*Fem.Dim,Fem.Mesh.NNode*Fem.Dim);
Fnc = sparse(Fem.Mesh.NNode*Fem.Dim,1);  
Ftc = sparse(Fem.Mesh.NNode*Fem.Dim,1);  
SDF = Fem.system.Contact{1};

Y0 = Fem.Mesh.Node;
[U]  = meshfield(Fem,Fem.solver.sol.x);
[dU] = meshfield(Fem,Fem.solver.sol.dx);

% Y(:,1) = Y0(:,1) + U(:,1);
% Y(:,2) = Y0(:,2) + U(:,2);

Y  = Y0 + U;

eps = 1e-6;
d = SDF(Y);
Intersect = find((d(:,end))<eps);

if ~isempty(Intersect)
    I = Intersect;
    %id = I(1:size(Y(I,:),1),1);
    n1 = (SDF(Y(I,:)+repmat([eps,0],size(Y(I,:),1),1))-d(I,end))/eps;
    n2 = (SDF(Y(I,:)+repmat([0,eps],size(Y(I,:),1),1))-d(I,end))/eps;

    Ux = (d(I,end)).*n1(:,end);
    Uy = (d(I,end)).*n2(:,end);

    % normal vector, tangent vector
    N(:,1) =  n1(:,end);
    N(:,2) =  n2(:,end);
    T(:,2) =  n1(:,end);
    T(:,1) = -n2(:,end);

    V0 = dU(I,:);
    %vUy = (dUy(I,end));

    %V0 = [vUx,vUy];
    Vn = V0./sqrt(sum((V0.^2),2) + 1e-2);

    Ffric = (dot(Vn.',T.').').*T;

    % ReactionForceX = -omega*Ux;
    % ReactionForceY = -omega*Uy;

    RF = omegaN * [-Ux, -Uy];
    RF = sqrt(sum((RF.^2),2));

    % modify forcing
    Ix = 2*I(1:size(Y(I,:),1),1)-1;
    Iy = 2*I(1:size(Y(I,:),1),1);

    % eDof = sort([Ix;Iy]);
    % I = repmat(eDof,1,numel(I)); 
    % J = I';

    Fnc(Ix,1) = -omegaN*Ux;
    Fnc(Iy,1) = -omegaN*Uy;
    Ftc(Ix,1) = -omegaT*(RF).*Ffric(:,1);% + 0*Dmod*V0(:,1);
    Ftc(Iy,1) = -omegaT*(RF).*Ffric(:,2);% + 0*Dmod*V0(:,2);
    
    index = 0;
    NDof  = Fem.Dim;
    I = zeros(numel(Ix)*NDof^2,1);
    J = zeros(numel(Ix)*NDof^2,1);
    K = zeros(numel(Ix)*NDof^2,1);
    T = zeros(numel(Ix)*NDof^2,1);

    for ii = 1:numel(Ix)
        KK = N(ii,:) .' * N(ii,:);
        TT = sign(RF(ii)) * T(ii,:) .' * T(ii,:);
        eDof = [Ix(ii);Iy(ii)];
        II = repmat(eDof,1,Fem.Dim);
        JJ = II';

        ind = index+1:index+NDof^2;
        I(ind) = II(:);
        J(ind) = JJ(:);
        K(ind) = KK(:);
        T(ind) = TT(:);

        index = index + NDof^2;
    end

    Knc = omegaN*sparse(I,J,K, ...
          Fem.Mesh.NNode*Fem.Dim,Fem.Mesh.NNode*Fem.Dim);

    Ktc = omegaT*sparse(I,J,T, ...
        Fem.Mesh.NNode*Fem.Dim,Fem.Mesh.NNode*Fem.Dim);

end

end

% % if isa(SDF,'Sdf')
% %     SDF = @
% % end

% %Move = Fem.Contact{2};
% % Emod = Fem.Material.getModulus();
% % Rmod = Fem.Material.getContactReaction();
% % Cmod = Fem.Material.getContactFriction();
% % Dmod = Fem.Material.getContactDamping();

% Y0 = Fem.Mesh.Node;
% [U]  = meshfield(Fem,Fem.solver.x);
% [dU] = meshfield(Fem,Fem.solver.dx);

% % Y(:,1) = Y0(:,1) + Ux - clamp(Fem.Time,0,1)*Move(1);
% % Y(:,2) = Y0(:,2) + Uy - clamp(Fem.Time,0,1)*Move(2);

% Y(:,1) = Y0(:,1) + U(:,1);
% Y(:,2) = Y0(:,2) + U(:,2);

% if Fem.Dim == 3
%     %Y(:,3) = Y0(:,3) + Uz - clamp(Fem.Time,0,1)*Move(3);
%     Y(:,2) = Y0(:,2) + U(:,3);
% end

% eps = 1e-6;
% d = SDF(Y);

% % if isa(SDF,'Sdf')
% %     d = SDF(Y);
% % else
% %     d = SDF.eval(Y);
% %     d = d(:,end);
% % end

% I = find((d(:,end))<eps);

% if ~isempty(I)
%         id = I(1:size(Y(I,:),1),1);
%         n1 = (SDF(Y(I,:)+repmat([eps,0],size(Y(I,:),1),1))-d(I,end))/eps;
%         n2 = (SDF(Y(I,:)+repmat([0,eps],size(Y(I,:),1),1))-d(I,end))/eps;

%         Ux = (d(I,end)).*n1(:,end);
%         Uy = (d(I,end)).*n2(:,end);

%         % normal vector, tangent vector
%         N(:,1) =  n1(:,end);
%         N(:,2) =  n2(:,end);
%         T(:,2) =  n1(:,end);
%         T(:,1) = -n2(:,end);

%         vUx = (dUx(I,end));
%         vUy = (dUy(I,end));

%         V0 = [vUx,vUy];
%         Vn = V0./sqrt(sum((V0.^2),2) + 1e-2);

%         Ffric = (dot(Vn.',T.').').*T;

%         %IPC = @(x,b) log(x/b).*(2*b - 2*x) - ((b - x).^2)./x;
%         %CUB = @(x) x.^3;
%         LIN = @(x) x;

%         Ux  = LIN(d(I,end).*n1(:,end));
%         Uy  = LIN(d(I,end).*n2(:,end));

%         vUx = LIN(d(I,end).*n1(:,end));
%         vUy = LIN(d(I,end).*n2(:,end));

%         ReactionForceX = -Rmod*Ux;
%         ReactionForceY = -Rmod*Uy;

%         RF = [ReactionForceX, ReactionForceY];
%         RF = sqrt(sum((RF.^2),2));

%         % modify forcing
%         Ix = 2*I(1:size(Y(I,:),1),1)-1;
%         Iy = 2*I(1:size(Y(I,:),1),1);

%         Fnc(Ix,1) = ReactionForceX;
%         Fnc(Iy,1) = ReactionForceY;

%         Ftc(2*id-1,1) = -Cmod*(RF).*Ffric(:,1) + Dmod*vUx;
%         Ftc(2*id,1)   = -Cmod*(RF).*Ffric(:,2) + Dmod*vUy;

% end
