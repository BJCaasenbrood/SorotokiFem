function Material = Ecoflex0030(D)
if nargin < 1, D = 2; end
Material = Yeoh('C1',17e-3,'C2',-0.2e-3,'C3',0.023e-3,...
    'D1',D,'D2',D,'D3',D);

% Material = YeohMaterial('C1',2.6776e-3 ,'C2',3.0716e-6,'C3',-3.9031e-9,...
%      'D1',D,'D2',D,'D3',D);

Material.params.Rho  = 1070e-12;
Material.params.Zeta = 4;
Material.contact.TangentFriction = 0.75;
Material.contact.NormalReaction = 5.75;
end

