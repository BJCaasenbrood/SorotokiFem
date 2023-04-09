function F = addBasicLoadsFem(Fem, F)

NLoad = size(Fem.system.Load,1);

for ii = 1:Fem.Dim
    F(Fem.Dim*Fem.system.Load(1:NLoad,1)+(ii-Fem.Dim),1) = ...
        Fem.system.Load(1:NLoad,1+ii);
end

end

