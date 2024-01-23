function Fem = addGravity(Fem,varargin)
    if isempty(varargin)
       if Fem.Dim == 2
           varargin{1} = [0,-9.81e3].';
       else
           varargin{1} = [0,0,-9.81e3].'; 
       end
    end
    Fem = addGravityFem(Fem,[],varargin{1:end});
end