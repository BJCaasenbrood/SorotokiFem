function Fem = addContact(Fem,varargin)
    if isa(varargin{1},'Sdf')
        sdf = varargin{1};
        varargin{1} = @(x) sdf.eval(x);
        B1 = box2node(sdf.BdBox); 
        B2 = box2node(Fem.BdBox);
        Fem.BdBox = boxhull([B1;B2],mean(abs(Fem.BdBox))/1e1);
        Fem.system.ContactSDF = sdf;
    else
        error('Input for Fem.addContact must be Sdf class');
    end
       
    if numel(varargin) < 2
        varargin{2} = [0,0]; 
    end
    
    Fem = addContactFem(Fem,varargin{1:end});

end