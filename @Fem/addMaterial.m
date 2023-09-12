% addMaterial - Add a material to the Fem structure.
%
%   Fem = addMaterial(Fem, Material) adds a material to the Fem structure.
%   The material is specified by the input argument Material.
%
%   Input:
%       - Fem: The Fem structure.
%       - Material: The material to be added.
%
%   Example:
%       Fem = addMaterial(Fem, NeoHookean)
%       Fem = addMaterial(Fem, Yeoh)

function Fem = addMaterial(Fem,Material)
    
    if isempty(Fem.materials.MatElem)
        Fem.materials.Material = {Material};
        Fem.materials.MatElem = ones(Fem.Mesh.NElem,1);
    else
        Fem.materials.Material{end+1} = Material;
    end
    Fem.materials.NMat = Fem.materials.NMat + 1;
end