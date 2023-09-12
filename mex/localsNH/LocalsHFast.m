function [Fe,Fb,Qe,Me,Ce,Ke,Kte,Svme,SS,EE,Te,Ve,Vge,Tke,Re,Ue] = ...
        LocalsHFast(eNode,eDof,dV,Rb,...
        Dim,...   % Fem.Dim
        Node0,... % Node0
        Ns,...     % shpfnc for Fem.ShapeFnc{nn}
        dNdxis,... % derive shpfnc for Fem.ShapeFnc{nn}
        W,...     % weights shapefnc
        Utmp,...  % current displacements
        dUtmp,... % current velocities
        Rho,...   % density
        Zeta,...  % dampings
        Grav,...  % gravity
        E0,...     
        Nu0...
        )
% get order
nn = length(eNode);
mm = Dim;

% get gauss points and weights
%W   = Fem.ShapeFnc{nn}.W;
Fe  = zeros(mm*nn,1);
Fb  = zeros(mm*nn,1);
Te  = zeros(mm*nn,1);
Me  = zeros(mm*nn,mm*nn);
Ce  = zeros(mm*nn,mm*nn);
Ke  = zeros(mm*nn,mm*nn);
Kte = zeros(mm*nn,mm*nn);
SGP = zeros(length(W),6);
EGP = zeros(length(W),6);
RRe = zeros(3);
UUe = zeros(3);
Qe  = ones(nn,1);
Ve  = 0;
Vge = 0;

Nshp = length(Ns(:,:,1));
NNe  = zeros(length(W),Nshp);

if mm == 2, Et = [dV;dV;0];
else, Et = [dV;dV;dV;0;0;0];
end

% get displacement field
Delta  = Utmp(eDof,:);
dDelta = dUtmp(eDof,:);

% F  = eye(3);
De = E0/(1-Nu0^2)*[1 Nu0 0;Nu0 1 0;0 0 (1-Nu0)/2];

% quadrature loop
for q = 1:length(W)
    
    % extract shape-functions
    N     = Ns(:,:,q);
    dNdxi = dNdxis(:,:,q);
    
    if numel(N) == 3 && Dim == 3
        T = cross(dNdxi(:,1),dNdxi(:,2));
        dNdxi = [dNdxi,(T)];
    end
    
    J0    = Node0(eNode,:).'*dNdxi;
    dNdx  = dNdxi/J0;
    dJ    = det(J0);
    
    % deformation gradient   
    F = DeformationGradient(Delta,dNdx,Dim);
    
    % get internal stress matrix
    % [S0, D0, Psi] = PiollaStress(Mu,Lambda,F);
    
    % voigt-notation vectorize
    % S = VoightNotation(S0);
    
    % reduced isotropic matrices
    S = De * Bnl * Delta;
    [~, De, Ge] = IsotropicReduction(D0,S,Dim);

    % polar decompostion
    % [Fiso, Fvol, ~] = PolarDecomposition(F);
    
    % nonlinear strain-displacement operator
    [Bnl,Bg,NN,tau] = NonlinearStrainOperatorFast(N,dNdx,F);
    
    % % local elemental rotation
    % RRe = RRe + Fiso/nn;
    % UUe = UUe + Fvol/nn;
    % NN

    % internal force vector
    Fe = Fe + tau*W(q)*Bnl.' * De * Bnl * Delta * dJ;
    
    % (graviational) body force vector
    Fb = Fb + tau*W(q)*Rho*(NN.')*Grav(:)*dJ;
    
    % lineararized stiffness matrix
    Kte = Kte + tau*W(q)*(Bnl.'*De*Bnl + Bg.'*Ge*Bg)*dJ;
    
    % mass matrix
    Me = Me + tau*W(q)*Rho*(NN.')*NN*dJ;
    
    % dampings matrix
    Ce = Ce + Zeta*Me;
    %Ce = Ce + Zeta*Kte;
    
    % thermal expansion force
    Te = Te + tau*W(q)*Bnl.'*De*Et*dJ;
    
    % elemental potential energy
    Ve = Ve  + tau*W(q)*Psi*dJ;
    
    % graviational energy
    Vgetmp = tau*W(q)*((N.'*Node0(eNode,:)) + ...
        (NN*Delta).')*Rho*Grav(:)*dJ;

    Vge = Vge - Vgetmp(1);
    
    % lagrangian strain
    Elagran = (1/2)*(F.'*F - eye(3));
    % Elagran = (1/2)*(zeros(3));
    
    % % true stress and lagrangian strain
    % SGP(q,:) = VoightNotation((1/det(F))*F*S0*(F.'));
    % EGP(q,:) = VoightNotation(Elagran);
    
    % construct shape functions
    NNe(((q-1)*Nshp + 1):(q*Nshp)) = N(:).';
end

% compute elemental kinetic energy
Tke = 0.5*(dDelta).'*Me*(dDelta);    

%compute elemental rotation matrix
[Ur,~,Vr] = svd(RRe);
Re = (Ur*Vr.');
Ue = UUe;

SS = NNe.'*SGP;
EE = NNe.'*EGP;
[Svm, ~] = VonMises(SS(:,1),SS(:,2),SS(:,3),...
                    SS(:,4),SS(:,5),SS(:,6));
Svme = Svm(:); 

end

%---------------------------------------------------------- polar decompose
function [R,S,V] = PolarDecomposition(F)
C = F.'*F;
[Q0, lambdasquare] = eig(C);

lambda = sqrt(diag((lambdasquare))); 
Uinv   = repmat(1./lambda',size(F,1),1).*Q0*Q0.';

R = real(F*Uinv);
S = real(R.'*F);
V = real(F*R.');
end

%------------------------------------------------ nonlinear strain operator
function [Bn,Bg,NN,tau] = NonlinearStrainOperatorFast(N,dNdx,F)
nn = length(N);
mm = size(dNdx,2);
zz = mm*nn;

NN = zeros(mm,zz);

id1 = 1:mm:zz;
id2 = 2:mm:zz;

NN(1,id1) = N.';
NN(2,id2) = N.';
dNdxX = dNdx(:,1).';
dNdxY = dNdx(:,2).';

if mm == 3
    NN(3,3:mm:zz) = N.';
    dNdxZ = dNdx(:,3).';
else
    dNdxZ = 0;
end

Bn = zeros((mm-1)*3,zz);
Bg = zeros((mm-1)*4+(mm-2),zz);

if mm == 2 % 2-dimensional
    Bn(1,id1) = dNdxX*F(1,1);
    Bn(1,id2) = dNdxX*F(2,1);
    Bn(2,id1) = dNdxY*F(1,2);
    Bn(2,id2) = dNdxY*F(2,2);
    Bn(3,id1) = dNdxX*F(1,2) + dNdxY*F(1,1);
    Bn(3,id2) = dNdxX*F(2,2) + dNdxY*F(2,1);
    
    Bg(1,id1) = dNdxX;
    Bg(2,id1) = dNdxY;
    Bg(3,id2) = dNdxX;
    Bg(4,id2) = dNdxY;
    
else % 3-dimensional
    Bn(1,1:mm:mm*nn) = dNdxX*F(1,1);
    Bn(1,2:mm:mm*nn) = dNdxX*F(2,1);
    Bn(1,3:mm:mm*nn) = dNdxX*F(3,1);
    Bn(2,1:mm:mm*nn) = dNdxY*F(1,2);
    Bn(2,2:mm:mm*nn) = dNdxY*F(2,2);
    Bn(2,3:mm:mm*nn) = dNdxY*F(3,2);
    Bn(3,1:mm:mm*nn) = dNdxZ*F(1,3);
    Bn(3,2:mm:mm*nn) = dNdxZ*F(2,3);
    Bn(3,3:mm:mm*nn) = dNdxZ*F(3,3);
    Bn(4,1:mm:mm*nn) = dNdxX*F(1,2) + dNdxY*F(1,1);
    Bn(4,2:mm:mm*nn) = dNdxX*F(2,2) + dNdxY*F(2,1);
    Bn(4,3:mm:mm*nn) = dNdxX*F(3,2) + dNdxY*F(3,1);
    Bn(5,1:mm:mm*nn) = dNdxY*F(1,3) + dNdxZ*F(1,2);
    Bn(5,2:mm:mm*nn) = dNdxY*F(2,3) + dNdxZ*F(2,2);
    Bn(5,3:mm:mm*nn) = dNdxY*F(3,3) + dNdxZ*F(3,2);
    Bn(6,1:mm:mm*nn) = dNdxX*F(1,3) + dNdxZ*F(1,1);
    Bn(6,2:mm:mm*nn) = dNdxX*F(2,3) + dNdxZ*F(2,1);
    Bn(6,3:mm:mm*nn) = dNdxX*F(3,3) + dNdxZ*F(3,1);
    
    Bg(1,1:mm:mm*nn) = dNdxX;
    Bg(2,1:mm:mm*nn) = dNdxY;
    Bg(3,1:mm:mm*nn) = dNdxZ;
    Bg(4,2:mm:mm*nn) = dNdxX;
    Bg(5,2:mm:mm*nn) = dNdxY;    
    Bg(6,2:mm:mm*nn) = dNdxZ;   
    Bg(7,3:mm:mm*nn) = dNdxX;
    Bg(8,3:mm:mm*nn) = dNdxY;    
    Bg(9,3:mm:mm*nn) = dNdxZ;     
end

tau = 1;

end 
%---------------------------------------------------------- polar decompose
function F = DeformationGradient(U,dNdx,Dim)
nn = round(length(U)/Dim);
UU = zeros(nn,Dim);
id1 = round(1:Dim:Dim*nn);
id2 = round(2:Dim:Dim*nn);
UU(:,1) = U(id1(:));
UU(:,2) = U(id2(:));

if Dim == 2
    F0 = (UU.') * dNdx;%(dNdx'*UU)';
    F = [F0(1,1)+1,F0(1,2),0; F0(2,1),F0(2,2)+1,0;0,0,1];
else
    id3 = round(3:Dim:Dim*nn).';
    UU(:,3) = U(id3);
    F = (dNdx'*UU)' + eye(3);
end
end
%------------------------------------------------ nonlinear strain operator
function [S, D, G] = IsotropicReduction(D0,S0,Dim)

if Dim == 2
    G = zeros(4,4);
    D   = [D0(1,1), D0(1,2),       0;
           D0(2,1), D0(2,2),       0;
                 0,       0, D0(4,4)];
    SIG = [S0(1), S0(4); S0(4), S0(2)];
    S   = [S0(1); S0(2); S0(4)]; 
    
    G(1:2,1:2) = SIG;
    G(3:4,3:4) = SIG;
else
    G = zeros(9,9);
    D = D0;
    S = S0;
    SIG = [S0(1), S0(4), S0(6);
           S0(4), S0(2), S0(5);
           S0(6), S0(5), S0(3)];
       
    G(1:3,1:3) = SIG;
    G(4:6,4:6) = SIG;
    G(7:9,7:9) = SIG;
end

end
%------------------------------------------------ nonlinear strain operator
function [Svm, Svmm] = VonMises(S11,S22,S33,S12,S23,S13)
s11 = S11; s22 = S22; s33 = S33; s12 = S12; s23 = S23; s13 = S13;
Svm = sqrt(0.5*((s11-s22).^2 + (s22-s33).^2 + (s33-s11).^2 ...
    + 6*(s12.^2 + s23.^2 + s13.^2)));
Svmm = mean(Svm); 
end
%--------------------------------------------------- Kelvin-voight notation
function Sv = VoightNotation(S)
Sv = [S(1,1); S(2,2); S(3,3); S(1,2); S(2,3); S(1,3)]; 
end