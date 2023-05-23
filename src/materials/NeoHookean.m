classdef NeoHookean

    properties (Access = public)
        Type;
        params;
        contact;

        % Mu; 
        % Lambda;
        % E;
        % Nu;
        % % Rho  = 970e-12;
        % % Zeta = 0.4;
        % % Rr   = 0.2;      
        % % Cfr  = 0.2;
        % Density                = 970e-12;
        % Damping                = 0.4;
        % ContactNormalDamping   = 0.05;
        % ContactNormalReaction  = 0.1;
        % ContactTangentFriction = 0.1;
    end
    
    properties (Access = private)
        Alpha;
    end

    
%--------------------------------------------------------------------------
methods  
%--------------------------------------------------------------- Mesh Class
function obj = NeoHookean(varargin) 
    
    obj.Type = 'NeoHookean';
    obj.params  = struct;
    obj.contact = struct;

    if isempty(varargin)
        E0  = 1;
        Nu0 = 0.45;
    elseif numel(varargin{1}) == 2
        C = varargin{1};
        E0 = C(1); Nu0 = C(2);
        if numel(varargin) > 1
            varargin = varargin(2:end);
        end
    else
        E0 = varargin{1};
        Nu0 = varargin{2};     
        if numel(varargin) > 2
            varargin = varargin(3:end);
        end   
    end
    
    obj.params.E      = E0;
    obj.params.Nu     = Nu0;
    obj.params.Lambda = (Nu0*E0)/((1+Nu0)*(1-2*Nu0));
    obj.params.Mu     = E0/(2*(1+Nu0));
    obj.params.Rho    = 970e-12;
    obj.params.Zeta   = 0.4;

    obj.contact.NormalDamping   = 0.05;
    obj.contact.NormalReaction  = 0.5;
    obj.contact.TangentFriction = 0.3;

    for ii = 1:2:length(varargin)
        if isfield(obj.contact,varargin{ii})
            obj.contact.(varargin{ii}) = varargin{ii+1};
        elseif isfield(obj.params,varargin{2*ii-1})
            obj.params.(varargin{ii}) = varargin{ii+1};
        end
    end

end
%---------------------------------------------------------------------- get     
function varargout = get(NeoHookean,varargin)
    if nargin > 1
        varargout{nargin-1,1} = [];
        for ii = 1:length(varargin)
            varargout{ii,1} = NeoHookean.(varargin{ii});
        end
    else
        varargout = NeoHookean.(varargin);
    end
end   
%---------------------------------------------------------------------- get     
function y = getModulus(NeoHookeanMaterial)
y = NeoHookeanMaterial.params.E;
end
%---------------------------------------------------------------------- get     
function y = getContactReaction(NeoHookean)
y = NeoHookean.contact.NormalReaction*...
    getModulus(NeoHookean);
end
%---------------------------------------------------------------------- get     
function y = getContactFriction(NeoHookean)
y = NeoHookean.contact.TangentFriction;
end
%---------------------------------------------------------------------- get     
function y = getContactDamping(NeoHookean)
y = NeoHookean.ContactNormalDamping...
    *getModulus(NeoHookean);
end
%---------------------------------------------------------------------- get    
function [S,S11,S33] = uniaxial(NeoHookeanMaterial,x)
S   = zeros(numel(x),1);
S11 = zeros(numel(x),1);
S33 = zeros(numel(x),1);
for ii = 1:numel(x)
    F   = diag([x(ii),sqrt(1/x(ii)),sqrt(1/x(ii))]);
    PK2 = PiollaStress(NeoHookeanMaterial,F);
    SIG = (1/det(F))*F*PK2*(F.');
    
    S(ii,1)   = SIG(1,1) - SIG(3,3);
    S11(ii,1) = SIG(1,1);
    S33(ii,1) = SIG(3,3);
end
end
%---------------------------------------------------------------------- set
function NeoHookeanMaterial = set(NeoHookeanMaterial,varargin)
    for ii = 1:2:length(varargin)
        NeoHookeanMaterial.(varargin{ii}) = varargin{ii+1};
    end
end
%------------------------------ 2nd Piolla stress tensor for Neo-Hookean
function [S, D, P] = PiollaStress(NeoHookeanMaterial,F,~)
%Se = 2nd PK stress [S11, S22, S33, S12, S23, S13];
C100 = NeoHookeanMaterial.params.Mu/2; 
K    = NeoHookeanMaterial.params.Lambda/2;

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

I3EE = [ 0     4*C3  4*C2  0    -4*C5  0;
         4*C3  0     4*C1  0     0    -4*C6;
         4*C2  4*C1  0    -4*C4  0     0;
         0     0    -4*C4 -2*C3  2*C6  2*C5;
        -4*C5  0     0     2*C6 -2*C1  2*C4;
         0    -4*C6  0     2*C5  2*C4 -2*C2];
%
W1 = X23*I3^(-X12);    W2 = X89*I1*I3^(-X43); W3 = X13*I1*I3^(-X43);
W8 = I3^(-X12);        W9 = X12*I3^(-X12);
%
J1EE = -W1*(J1E*J3E' + J3E*J1E') + W2*(J3E*J3E') - W3*I3EE;
J3EE = -W8*(J3E*J3E') + W9*I3EE;
%
D = C100*J1EE + K*(J3E*J3E') + K*J3M1*J3EE;

end
%---------------------------------------------------------------------- set
function E = Emod(NeoHookeanMaterial)
   E = NeoHookeanMaterial.E;
end

end
end

