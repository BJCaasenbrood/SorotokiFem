function obj = materialoptions
    obj = NamedTuple('Material',{});
    obj.add('NMat',0);
    obj.add('Material',{});
    obj.add('ElemMat',[]);
    obj.add('ColorMap',cmap_materials);
end

