function Fem = assembleDeformationGrad(Fem)

    list = 1:Fem.Mesh.NElem;
    F2V = Fem.Mesh.geometry.FaceToNode;

    % for eacth element
    for ii = 1:Fem.Mesh.NNode
        Rtmp = zeros(3);
        Qtmp = zeros(3);
        ElemFromNode = list(F2V(:, ii) > 0);
        % 
        if numel(ElemFromNode) > 1
        %     A = Fem.system.F.ElemStr(ElemFromNode);
        %     B = Fem.system.F.ElemRot(ElemFromNode);
        % 
        %     Qtmp = sum(reshape(cell2mat(A),size(A{1},1),[],numel(ElemFromNode)),3);
        %     Rtmp = sum(reshape(cell2mat(B),size(B{1},1),[],numel(ElemFromNode)),3);
        % 
             Qtmp = sum(cellfun(@(x) sum(x(:)),Fem.system.F.ElemStr(ElemFromNode)));
             Rtmp = sum(cellfun(@(x) sum(x(:)),Fem.system.F.ElemRot(ElemFromNode)));
         else
             Qtmp = Fem.system.F.ElemStr{ElemFromNode};
             Rtmp = Fem.system.F.ElemRot{ElemFromNode};
         end

        %Qtmp = cellfun(@sum,Fem.system.F.ElemStr{ElemFromNode});
        %Rtmp = cellfun(@sum,Fem.system.F.ElemRot{ElemFromNode});
        % 
        % for jj = ElemFromNode
        %     Qtmp = Qtmp + Fem.system.F.ElemStr{jj};
        %     Rtmp = Rtmp + Fem.system.F.ElemRot{jj};
        % end

        [Ur, ~, Vr] = svd(Rtmp);
        R{ii} = (Ur * Vr.');
        Q{ii} = Qtmp / numel(ElemFromNode);
    end

    Fem.system.Rotation = R;
    Fem.system.Stretch  = Q;
    
end
