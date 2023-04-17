function Fem = assembleDeformationGrad(Fem)

    list = 1:Fem.Mesh.NElem;
    F2V = Fem.Mesh.geometry.FaceToNode;

    % for eacth element
    for ii = 1:Fem.Mesh.NNode
        Rtmp = 0;
        Qtmp = 0;
        ElemFromNode = list(F2V(:, ii) > 0);

        for jj = ElemFromNode
            Qtmp = Qtmp + Fem.system.F.ElemStr{jj};
            Rtmp = Rtmp + Fem.system.F.ElemRot{jj};
        end

        [Ur, ~, Vr] = svd(Rtmp);
        R{ii} = (Ur * Vr.');
        Q{ii} = Qtmp / numel(ElemFromNode);
    end

    Fem.system.Rotation = R;
    Fem.system.Stretch  = Q;
    
end
