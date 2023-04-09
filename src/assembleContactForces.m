function F = assembleContactForces(Fem, F)

    SDF  = Fem.Contact{1};
    Move = Fem.Contact{2};
    Emod = Fem.Material.getModulus();
    Rmod = Fem.Material.getContactReaction();
    Cmod = Fem.Material.getContactFriction();
    Dmod = Fem.Material.getContactDamping();

    Y0 = Fem.Mesh.Node0;
    [Ux,Uy,~,Uz]    = DisplacementField(Fem,Fem.Utmp);
    [dUx,dUy,~,dUz] = DisplacementField(Fem,Fem.dUtmp);
    Y(:,1) = Y0(:,1) + Ux - clamp(Fem.Time,0,1)*Move(1);
    Y(:,2) = Y0(:,2) + Uy - clamp(Fem.Time,0,1)*Move(2);

    if Fem.Dim == 3
        Y(:,3) = Y0(:,3) + Uz - clamp(Fem.Time,0,1)*Move(3);
    end

    eps = 1e-6;

    d = SDF(Y);
    I = find((d(:,end))<eps);

    if ~isempty(I)
        if Fem.Dim == 3
            n1 = (SDF(Y(I,:)+repmat([eps,0,0],size(Y(I,:),1),1))-d(I,end))/eps;
            n2 = (SDF(Y(I,:)+repmat([0,eps,0],size(Y(I,:),1),1))-d(I,end))/eps;
            n3 = (SDF(Y(I,:)+repmat([0,0,eps],size(Y(I,:),1),1))-d(I,end))/eps;

            Ux = (d(I,end)).*n1(:,end);
            Uy = (d(I,end)).*n2(:,end);
            Uz = (d(I,end)).*n3(:,end);

            Fc(3*I(1:size(Y(I,:),1),1)-2,1) = -Rmod*(Ux).^1;
            Fc(3*I(1:size(Y(I,:),1),1)-1,1) = -Rmod*(Uy).^1;
            Fc(3*I(1:size(Y(I,:),1),1),1)   = -Rmod*(Uz).^1;

            C(3*I(1:size(Y(I,:),1),1)-1,2*I(1:size(Y(I,:),1),1)-1) = ...
                C(2*I(1:size(Y(I,:),1),1)-1,2*I(1:size(Y(I,:),1),1)-1) + Cmod;

            C(3*I(1:size(Y(I,:),1),1),2*I(1:size(Y(I,:),1),1)) = ...
                C(2*I(1:size(Y(I,:),1),1),2*I(1:size(Y(I,:),1),1)) + Cmod;

            C(3*I(1:size(Y(I,:),1),1),2*I(1:size(Y(I,:),1),1)) = ...
                C(2*I(1:size(Y(I,:),1),1),2*I(1:size(Y(I,:),1),1)) + Cmod;

        else
            id = I(1:size(Y(I,:),1),1);
            n1 = (SDF(Y(I,:)+repmat([eps,0],size(Y(I,:),1),1))-d(I,end))/eps;
            n2 = (SDF(Y(I,:)+repmat([0,eps],size(Y(I,:),1),1))-d(I,end))/eps;

            Ux = (d(I,end)).*n1(:,end);
            Uy = (d(I,end)).*n2(:,end);

            % normal vector, tangent vector
            N(:,1) =  n1(:,end);
            N(:,2) =  n2(:,end);
            T(:,2) =  n1(:,end);
            T(:,1) = -n2(:,end);

            vUx = (dUx(I,end));
            vUy = (dUy(I,end));

            V0 = [vUx,vUy];
            Vn = V0./sqrt(sum((V0.^2),2) + 1e-2);

            Ffric = (dot(Vn.',T.').').*T;

            %IPC = @(x,b) log(x/b).*(2*b - 2*x) - ((b - x).^2)./x;
            %CUB = @(x) x.^3;
            LIN = @(x) x;

            Ux  = LIN(d(I,end).*n1(:,end));
            Uy  = LIN(d(I,end).*n2(:,end));

            vUx = LIN(d(I,end).*n1(:,end));
            vUy = LIN(d(I,end).*n2(:,end));

            ReactionForceX = -Rmod*Ux;
            ReactionForceY = -Rmod*Uy;

            RF = [ReactionForceX, ReactionForceY];
            RF = sqrt(sum((RF.^2),2));

            % modify forcing
            Ix = 2*I(1:size(Y(I,:),1),1)-1;
            Iy = 2*I(1:size(Y(I,:),1),1);

            Fnc(Ix,1) = ReactionForceX;
            Fnc(Iy,1) = ReactionForceY;

            Ftc(2*id-1,1) = -Cmod*(RF).*Ffric(:,1) + Dmod*vUx;
            Ftc(2*id,1)   = -Cmod*(RF).*Ffric(:,2) + Dmod*vUy;

        end
    end
end

