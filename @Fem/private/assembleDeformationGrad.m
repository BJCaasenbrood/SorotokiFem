function Fem = assembleDeformationGrad(Fem)

    list = 1:Fem.Mesh.NElem;
    F2V = Fem.Mesh.geometry.FaceToNode;

    % for eacth element
    R = cell(Fem.Mesh.NNode, 1);
    Q = cell(Fem.Mesh.NNode, 1);

    for ii = 1:Fem.Mesh.NNode

        ElemFromNode = list(F2V(:, ii) > 0);

        if numel(ElemFromNode) > 1
             Qtmp = sum(cellfun(@(x) sum(x(:)),Fem.system.F.ElemStr(ElemFromNode)));
             Rtmp = sum(cellfun(@(x) sum(x(:)),Fem.system.F.ElemRot(ElemFromNode)));
         else
             Qtmp = Fem.system.F.ElemStr{ElemFromNode};
             Rtmp = Fem.system.F.ElemRot{ElemFromNode};
         end

        [Ur, ~, Vr] = svd(Rtmp);
        R{ii} = (Ur * Vr.');
        Q{ii} = Qtmp / numel(ElemFromNode);
    end

    Fem.system.Rotation  = R;
    Fem.system.Stretch   = Q;
    
end
