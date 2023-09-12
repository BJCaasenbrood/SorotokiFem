function P = generateFilter(Fem,PS1)

    PS2 = Fem.Center0;
    PS3 = Fem.Node0;
    %ShapeFnc = Shapes.Fem.ShapeFnc;
    
    d = cell(size(PS1,1),1);    
    
    % get closest neighbors with element centers
    I = knnsearch(PS1,PS2);
    
    for ii = 1:size(PS1,1)
       XI = zeros(1,Fem.Dim);
       el = Fem.Element{I(ii)};
       Vsp = PS3(el,:);  % points spanned by element
       Vin = PS1(ii,:);  % point in element
     
       xi = Fem.ShapeFnc{numel(el)}.Xi;
    
       if Fem.Dim == 2
           XI(1) = griddata(Vsp(:,1),Vsp(:,2),xi(:,1),Vin(1),Vin(2),'natural');
           XI(2) = griddata(Vsp(:,1),Vsp(:,2),xi(:,2),Vin(1),Vin(2),'natural');
           
           if isnan(XI(1)) || isnan(XI(2))
               E = Vsp - Vin;
               [~,id] = min(sum(E.^2,2));
               XI(1) = griddata(Vsp(:,1),Vsp(:,2),xi(:,1),...
                   Vsp(id,1),Vsp(id,2),'natural');
               XI(2) = griddata(Vsp(:,1),Vsp(:,2),xi(:,2),...
                   Vsp(id,1),Vsp(id,2),'natural');
           end
           
       else
           
           XI(1) = griddata(Vsp(:,1),Vsp(:,2),Vsp(:,3),xi(:,1),...
               Vin(1),Vin(2),Vin(3),'natural');
           XI(2) = griddata(Vsp(:,1),Vsp(:,2),Vsp(:,3),xi(:,2),...
               Vin(1),Vin(2),Vin(3),'natural');
           XI(3) = griddata(Vsp(:,1),Vsp(:,2),Vsp(:,3),xi(:,3),...
               Vin(1),Vin(2),Vin(3),'natural');
           
           if isnan(XI(1)) || isnan(XI(2)) || isnan(XI(2))
               E = Vsp - Vin;
               [~,id] = min(sum(E.^2,2));
               XI(1) = griddata(Vsp(:,1),Vsp(:,2),Vsp(:,3),xi(:,1),...
                   Vsp(id,1),Vsp(id,2),Vsp(id,3),'natural');
               XI(2) = griddata(Vsp(:,1),Vsp(:,2),Vsp(:,3),xi(:,2),...
                   Vsp(id,1),Vsp(id,2),Vsp(id,3),'natural');
               XI(3) = griddata(Vsp(:,1),Vsp(:,2),Vsp(:,3),xi(:,3),...
                   Vsp(id,1),Vsp(id,2),Vsp(id,3),'natural');
           end
           
       end
       
       N = Fem.ShapeFnc{numel(el)}.fnc(XI);
       
       % assemble distance filter matrix based on ShpFnc N(s)
       if numel(el) == 3
        d{ii} = [repmat(ii,numel(el),1),[el(2);el(1);el(3)],(N)];
       else
        d{ii} = [repmat(ii,numel(el),1),[(el(:))],(N)];   
       end
       
    end
    
    d = cell2mat(d); 
    P = sparse(d(:,1),d(:,2),d(:,3),size(PS1,1),Fem.NNode);
    P = spdiags(1./sum(P,2),0,size(P,1),size(P,1))*P;
     
end