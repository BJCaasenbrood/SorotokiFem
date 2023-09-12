classdef YeohIncompressible

    properties (Access = public)
        Type;
        params;
        contact;
    end
    
    properties (Access = private)
        ID;
        SET;
        WGT;
        TensorCalc = false;
    end
   
%--------------------------------------------------------------------------
methods  
%--------------------------------------------------------------- Mesh Class
function obj = YeohIncompressible(varargin) 

    obj.Type   = 'Yeoh';
    obj.params = struct;
    obj.contact = struct; 
    
    [obj.ID, obj.SET, obj.WGT] = Tensor4IdSymmetric;

    obj.params.Rho  = 970e-12;
    obj.params.Zeta = 0.4;

    obj.params.C1 = 1;
    obj.params.C2 = 0;
    obj.params.C3 = 0;

    obj.params.D1 = 0;
    obj.params.D2 = 0;
    obj.params.D3 = 0;

    obj.contact.NormalDamping   = 0.05;
    obj.contact.NormalReaction  = 0.1;
    obj.contact.TangentFriction = 0.1;
    
    if ~isempty(varargin)
        if numel(varargin{1}) == 3
            C = varargin{1};
            obj.params.C1 = C(1); 
            obj.params.C2 = C(2); 
            obj.params.C3 = C(3);
        else
            for ii = 1:2:length(varargin)
                obj.params.(varargin{ii}) = varargin{ii+1};
            end
        end
    end


    
end
%---------------------------------------------------------------------- get     
function varargout = get(Yeoh,varargin)
    if nargin > 1
        varargout{nargin-1,1} = [];
        for ii = 1:length(varargin)
            varargout{ii,1} = Yeoh.(varargin{ii});
        end
    else
        varargout = Yeoh.(varargin);
    end
end       
%---------------------------------------------------------------------- set
function Yeoh = set(Yeoh,varargin)
    for ii = 1:2:length(varargin)
        Yeoh.(varargin{ii}) = varargin{ii+1};
    end
end  
%---------------------------------------------------------------------- set
function y = getModulus(Yeoh)
    y = 6*Yeoh.C1; 
end
%---------------------------------------------------------------------- get     
function y = getContactReaction(Yeoh)
y = Yeoh.contact.NormalReaction*...
    getModulus(Yeoh);
end
%---------------------------------------------------------------------- get     
function y = getContactFriction(Yeoh)
y = Yeoh.contact.TangentFriction;
end
%---------------------------------------------------------------------- get     
function y = getContactDamping(Yeoh)
y = Yeoh.contact.NormalDamping...
    *getModulus(Yeoh);
end
%---------------------------------------------------------------------- set
function [S,S11,S33,P] = uniaxial(Yeoh,x)
S   = zeros(numel(x),1);
P   = zeros(numel(x),1);
S11 = zeros(numel(x),1);
S33 = zeros(numel(x),1);
for ii = 1:numel(x)
    F   = diag([x(ii),sqrt(1/x(ii)),sqrt(1/x(ii))]);
    [PK2,~,psi] = PiollaStress(Yeoh,F);
    SIG = (1/det(F))*F*PK2*(F.');
    
    P(ii,1)   = psi;
    S(ii,1)   = SIG(1,1) - SIG(3,3);
    S11(ii,1) = SIG(1,1);
    S33(ii,1) = SIG(3,3);
end
end
%------------------------------ 2ND PIOLLA STRESSAND STIFFNESS FOR YEOH
function [S, D, P] = PiollaStress(Yeoh,F,~)
%Se = 2nd PK stress [S11, S22, S33, S12, S23, S13];
P = 0;
S = zeros(3,3);
C = F.'*F;
J = det(F);

YeohC = [Yeoh.params.C1,Yeoh.params.C2,Yeoh.params.C3];
I     = eye(3,3);
C11   = C(1,1); 
C22   = C(2,2); 
C33   = C(3,3);
I1    = C11 + C22 + C33;
%I1iso = I1;

for ii = 1:3
    P = P + YeohC(ii)*(I1 - 3)^(ii);
end

for ii = 1:3
    S = S + 2*(ii*YeohC(ii)*(I1 - 3)^(ii-1))*I;
end

alpha = 0;

for ii = 2:3
    alpha = alpha + ii*(ii-1)*YeohC(ii)*(I1 - 3)^(ii-2);
end

II3 = I;
TOa  = TensorOperation(II3,II3,true);

D = (4*alpha)*TOa;

end
%------------------------------ 2ND PIOLLA STRESSAND STIFFNESS FOR YEOH
function y = dWdI(YeohMaterial,I1)
c1 = YeohMaterial.C1; 
c2 = YeohMaterial.C2; 
c3 = YeohMaterial.C3;
y = c1 + 2*c2*(I1 -3) + 3*c3*(I1 -3).^2;
end

end

methods (Access = private)
 
end
end

function T = TensorOperation(A,B,Arg)
%#codegen

a11 = A(1,1); a12 = A(1,2); a13 = A(1,3);  
a21 = A(2,1); a22 = A(2,2); a23 = A(2,3);  
a31 = A(3,1); a32 = A(3,2); a33 = A(3,3);  

b11 = B(1,1); b12 = B(1,2); b13 = B(1,3);  
b21 = B(2,1); b22 = B(2,2); b23 = B(2,3);  
b31 = B(3,1); b32 = B(3,2); b33 = B(3,3);  

if Arg 
    T = [   a11*b11,   a22*b11,   a33*b11, 2*a12*b11, 2*a23*b11, 2*a13*b11;
   a11*b22,   a22*b22,   a33*b22, 2*a12*b22, 2*a23*b22, 2*a13*b22;
   a11*b33,   a22*b33,   a33*b33, 2*a12*b33, 2*a23*b33, 2*a13*b33;
 2*a11*b12, 2*a22*b12, 2*a33*b12, 4*a12*b12, 4*a23*b12, 4*a13*b12;
 2*a11*b23, 2*a22*b23, 2*a33*b23, 4*a12*b23, 4*a23*b23, 4*a13*b23;
 2*a11*b13, 2*a22*b13, 2*a33*b13, 4*a12*b13, 4*a23*b13, 4*a13*b13];
else
    T = [      a11*b11,           a21*b21,           a31*b31,             2*a11*b21,             2*a21*b31,             2*a11*b31;
           a12*b12,           a22*b22,           a32*b32,             2*a12*b22,             2*a22*b32,             2*a12*b32;
           a13*b13,           a23*b23,           a33*b33,             2*a13*b23,             2*a23*b33,             2*a13*b33;
 a11*b12 + a12*b11, a21*b22 + a22*b21, a31*b32 + a32*b31, 2*a11*b22 + 2*a12*b21, 2*a21*b32 + 2*a22*b31, 2*a11*b32 + 2*a12*b31;
 a12*b13 + a13*b12, a22*b23 + a23*b22, a32*b33 + a33*b32, 2*a12*b23 + 2*a13*b22, 2*a22*b33 + 2*a23*b32, 2*a12*b33 + 2*a13*b32;
 a11*b13 + a13*b11, a21*b23 + a23*b21, a31*b33 + a33*b31, 2*a11*b23 + 2*a13*b21, 2*a21*b33 + 2*a23*b31, 2*a11*b33 + 2*a13*b31];
end

end

function [id, set, W] = Tensor4IdSymmetric

id = [1,1; 2,2; 3,3; 1,2; 2,3; 1,3];

set = {[1,1], [1,2], [1,3], [1,4], [1,5], [1,6],...
       [2,1], [2,2], [2,3], [2,4], [2,5], [2,6],...
       [3,1], [3,2], [3,3], [3,4], [3,5], [3,6],...
       [4,1], [4,2], [4,3], [4,4], [4,5], [4,6],...
       [5,1], [5,2], [5,3], [5,4], [5,5], [5,6],...
       [6,1], [6,2], [6,3], [6,4], [6,5], [6,6]};
   
set = set(:);

W = kron([1,sqrt(2);sqrt(2),2],ones(3));%%
%W = kron([1,2;2,4],ones(3));
% W = [1,1,1,2,2,2;1,1,1,2,2,2;1,1,1,2,2,2;...
%      2,2,2,4,4,4;2,2,2,4,4,4;2,2,2,4,4,4];
end


