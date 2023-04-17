function [Fe] = ...
        LocalsNHFastElastic(eNode,eDof,dV,Rb,...
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
        Mu,...     
        Lambda...
        )
% get order
nn = length(eNode);
mm = Dim;

% get gauss points and weights
%W   = Fem.ShapeFnc{nn}.W;
Fe  = zeros(mm*nn,1);

if mm == 2, Et = [dV;dV;0];
else, Et = [dV;dV;dV;0;0;0];
end

% get displacement field
Delta  = Utmp(eDof,:);
dDelta = dUtmp(eDof,:);

% quadrature loop
for q = 1:length(W)
    
    % extract shape-functions
    N     = Ns(:,:,q);
    dNdxi = dNdxis(:,:,q);
    
    J0    = Node0(eNode,:).'*dNdxi;
    dNdx  = dNdxi/J0;
    dJ    = det(J0);
    
    % deformation gradient   
    F = DeformationGradient(Delta,dNdx,Dim);
    
    % get internal stress matrix
    [S0] = PiollaStress(Mu,Lambda,F);
    
    % voigt-notation vectorize
    S = VoightNotation(S0);
    
    % reduced isotropic matrices
    [Se] = IsotropicReduction(S,Dim);
    
    % nonlinear strain-displacement operator
    [Bnl] = NonlinearStrainOperatorFast(N,dNdx,F);
    
    % internal force vector
    Fe = Fe + W(q)*Bnl.'*Se*dJ;
end

end

function [S] = PiollaStress(Mu,Lambda,F)
%Se = 2nd PK stress [S11, S22, S33, S12, S23, S13];
C100 = Mu/2; 
K    = Lambda/2;

X12 = 1/2; 
X13 = 1/3; 
X23 = 2/3; 
X43 = 4/3; 
X89 = 8/9;
C = F.'*F;

C1=C(1,1); C2=C(2,2); C3=C(3,3); C4=C(1,2); C5=C(2,3); C6=C(1,3);

I1 = C1+C2+C3;
I3 = det(C);
J1 = I1*I3^(-X13);
J3 = sqrt(I3);
J3M1 = J3 - 1;
%
I1E = 2*[1,1,1,0,0,0]';
I3E = 2*[C2*C3-C5^2,  C3*C1-C6^2,  C1*C2-C4^2, ...
    C5*C6-C3*C4, C6*C4-C1*C5, C4*C5-C2*C6]';
%
W1 = I3^(-X13); W2 = X13*I1*I3^(-X43); W5 = X12*I3^(-X12);
%
J1E = W1*I1E - W2*I3E;
J3E = W5*I3E;
%

P = C100*(J1-3) + 0.5*K*(J3 - 1)^2;

Se = C100*J1E + K*J3M1*J3E;

S = [Se(1), Se(4), Se(6); 
     Se(4), Se(2), Se(5); 
     Se(6), Se(5), Se(3)];
end
%---------------------------------------------------------- polar decompose
function F = DeformationGradient(U,dNdx,Dim)
nn = round(length(U)/Dim);
UU = zeros(nn,Dim);
id1 = round(1:Dim:Dim*nn).';
id2 = round(2:Dim:Dim*nn).';
UU(:,1) = U(id1);
UU(:,2) = U(id2);

if Dim == 2
    F0 = (UU.')*dNdx;
    F = [F0(1,1)+1,F0(1,2),0; F0(2,1),F0(2,2)+1,0;0,0,1];
else
    id3 = round(3:Dim:Dim*nn).';
    UU(:,3) = U(id3);
    F = (dNdx'*UU)' + eye(3);
end
end

%------------------------------------------------ nonlinear strain operator
function [Bn] = NonlinearStrainOperatorFast(N,dNdx,F)
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

if mm == 2 % 2-dimensional
    Bn(1,id1) = dNdxX*F(1,1);
    Bn(1,id2) = dNdxX*F(2,1);
    Bn(2,id1) = dNdxY*F(1,2);
    Bn(2,id2) = dNdxY*F(2,2);
    Bn(3,id1) = dNdxX*F(1,2) + dNdxY*F(1,1);
    Bn(3,id2) = dNdxX*F(2,2) + dNdxY*F(2,1);
    
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
      
end

end 
%------------------------------------------------ nonlinear strain operator
function [S] = IsotropicReduction(S0,Dim)

if Dim == 2
    S   = [S0(1); S0(2); S0(4)]; 
else
    S = S0;
end

end
%--------------------------------------------------- Kelvin-voight notation
function Sv = VoightNotation(S)
Sv = [S(1,1); S(2,2); S(3,3); S(1,2); S(2,3); S(1,3)]; 
end