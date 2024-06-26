function msh = exportMesh(Fem,varargin)
    
    ISO = varargin{1};
    [~,I,UxUy] = showISO(Fem,ISO);
      
    B = UxUy;
    Xscale = (B(2)-B(1))/size(I.CData,2);
    Yscale = (B(4)-B(3))/size(I.CData,1);
    
    simplify_tol = varargin{2};
    
    img = I.CData >= 25;
    img = fliplr(img.');
    
    bnd = bwboundaries(img);
    
    c_cell0 = {};
    c_cell = {};
    
    for ii=1:length(bnd)
        bnd_tmp = bnd{ii};
        assert(all(bnd_tmp(1,:)==bnd_tmp(end,:)),'contour is not closed');
        c_cell0{ii} = bnd_tmp;
    end
    
    for ii=1:length(c_cell0)
        c_tmp = c_cell0{ii};
        c_red = decimatePoly(c_tmp,[simplify_tol, 2],false);
        if (nnz(c_red(:,1))>0)&&(nnz(c_red(:,2))>0)
            c_cell{end+1,1} = [Xscale*c_red(:,1), (Yscale)*c_red(:,2)];
        end
    end
    
    % create the 2d triangulation
    H = varargin{3};
    Tesselation = triangulationCreate(c_cell, H(1), H(2), H(3),'linear');
    
    msh = Mesh(Tesselation.Nodes.',Tesselation.Elements.');
    msh = msh.generate();
end