function h = showContactFem(Fem)

    sdf    = Fem.system.Contact{1};
    BD     = Fem.BdBox;

    [px,py] = meshgrid(linspace(BD(1),BD(2),80),...
                       linspace(BD(3),BD(4),80));
    Y = [px(:), py(:)];

    d = sdf(Y);
    D = reshape(d(:,end),[80,80]);

    color = [0.95,0.95,0.957];
    [~,h] = contourf(px,py,-D,[1e-3 1e-3], ...
        'EdgeColor',0.45*color,...
        'FaceColor',0.85*color,...
        'LineW',1.5);
end