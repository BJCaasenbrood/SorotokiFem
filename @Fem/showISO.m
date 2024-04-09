function [Fem,I,UxUy] = showISO(Fem, varargin)

    [V,X,Y,~] = former(Fem);

    % V = Fem.Topology;
    % X = Fem.TopologyGridX; 
    % Y = Fem.TopologyGridY; 

    [~,~,a] = size(V);

    if nargin < 3
        depth = floor(a/2);
    else
        depth = a * varargin{2}; 
    end

    X    = X(:,:,depth);
    Y    = Y(:,:,depth);
    SDF0 = V(:,:,depth);

    SDF = generateIsoCell(Fem,SDF0);

    if ~isempty(Fem.topology.CellRepetion)
        Rpt = Fem.topology.CellRepetion;
        SDF = repmat(SDF,Rpt(2),Rpt(1));
    else
        Rpt = [1, 1];
    end

    if ~isempty(Fem.topology.ReflectionPlane)
        Rp = Fem.topology.ReflectionPlane;
        
        if Rp(1) == 1
            Xf = flip(X) + (max(max(max(X))) - min(min(min(X))));
            X  = vertcat(Xf,X);
        end
        
        if Rp(1) == -1
            Xf = flip(X) - (max(max(max(X))) - min(min(min(X))));
            X  = horzcat(X,Xf);
        end
        
        if Rp(2) == -1
            Yf = flip(Y) - (max(max(max(Y))) - min(min(min(Y))));
            Y  = horzcat(Y,Yf);
        end
        
        if Rp(2) == 1
            Yf = flip(Y) + (max(max(max(Y))) - min(min(min(Y))));
            Y  = horzcat(Yf,Y);
        end
    end

    scaleX = 1; scaleY = 1;
    if ~isempty(Fem.Repeat)
        SDF0     = SDF;
        Instruct = Fem.Repeat;
        
        for ii = 1:length(Instruct)
            
            if Instruct(ii) == 1
                scaleX = scaleX + 1;
                SDF    = cat(2,SDF,SDF0);
            end
            
            if Instruct(ii) == 2
                scaleY = scaleY*2;
                SDF    = cat(1,SDF,SDF);
            end
        end
    end

    cla;
    Uxx  = scaleX*Rpt(1)*[min(min(min(X))) max(max(max(X)))];
    Uyy  = scaleY*Rpt(2)*[min(min(min(Y))) max(max(max(Y)))];
    UxUy = [0, Uxx(2)-Uxx(1), 0, Uyy(2)-Uyy(1)];

    I = smooth2a((SDF >= varargin{1})*255,3);

    I = image(rescale(Uxx),...
        ((max(Uyy)-min(Uyy))/(max(Uxx) - min(Uxx)))*rescale(Uyy),I);

    axis equal; axis off; 
    colormap(Fem.topology.ColorMap); 
    clim([0 1]);

    % background('w');
end

%
function [SDF,X,Y,Z] = former(Fem)
Res = 300;
Layers = 4;
Thickness = 0.2*(Fem.BdBox(2)-Fem.BdBox(1));

verts = Fem.Mesh.Node; 
fil   = full(Fem.topology.SpatialFilter);
V     = full(fem.Mesh.geometry.NodeToFace)*fil*fem.topology.sol.x;

% if Fem.VolumetricPressure
%     id = FindElements(Fem,'FloodFill',Fem,Fem.Density);
%     Pc = Fem.Mesh.Center;
%     S = zeros(Fem.NElem,1);
%     d = distancepointset(Pc, Pc(id,:), 0.5*Fem.FilterRadius );
%     ide = unique(d(:,2));
    
%     rhotmp =  fil*Fem.Density;
%     S(ide,1) = (rhotmp(ide)>Fem.VoidTolerance)*0.5;
%     Fill = zeros(Fem.NElem,1); 
%     Fill(id) = 0.5;
    
%     E = Fem.Mesh.get('NodeToFace')*...
%         fil*(S + Fem.Density*0.01);
%     V2 = clamp(Fem.Mesh.get('NodeToFace')*...
%         fil*(Fill + Fem.Density-0.1),0,1);
% end

% if isempty(Fem.Crop)
%     B = Fem.BdBox; 
% else
%     B = Fem.Crop;
% end

x = verts(:,1); 
y = verts(:,2);
dx = 0.01*(B(2) - B(1));
dy = 0.01*(B(4) - B(3));
xq = linspace(B(1)+dx,B(2)-dx,Res); 
yq = linspace(B(3)+dy,B(4)-dy,Res);

[xxq, yyq] = meshgrid(xq,yq);
P = griddata(x,y,double(V),xxq,yyq);
vrt = [xxq(:),yyq(:)];
Dist = Fem.Mesh.SDF(vrt); 
Dist = Dist(:,end);

tol = 0.1*sqrt((B(2)-B(1))*(B(4)-B(3))/length(vrt));

P = P(:); 
P(Dist > tol) = 0; 
P = (reshape(P,Res,Res)); 
V = cat(3,xxq,yyq,P);

% if Fem.VolumetricPressure
%     P = griddata(x,y,double(V2),xxq,yyq);
%     Dist = Fem.Mesh.SDF([xxq(:),yyq(:)]); 
%     Dist = Dist(:,end);
    
%     P = P(:); 
%     P(Dist > tol) = 0; 
%     P(isnan(P))=0;
%     P = (reshape(P,Res,Res));
    
%     GapFill = cat(3,xxq,yyq,P);
    
%     P = griddata(x,y,double(E),xxq,yyq);
%     Dist = Fem.Mesh.SDF([xxq(:),yyq(:)]); 
%     Dist = Dist(:,end);
    
%     P = P(:);  
%     P(isnan(P))=0;
%     P = (reshape(P,Res,Res));
    
%     EdgeFill = cat(3,xxq,yyq,P);
% end

SDF0 = V(:,:,3);
SDF0 = smooth2a(SDF0,5);
SDF  = repmat(SDF0,[1 1 Layers]);

% if Fem.VolumetricPressure
%     V       = GapFill;
%     SDF2    = V(:,:,3);
%     SDF2    = GaussianFilter(SDF2,10);
%     ZFiller = repmat(SDF2,[1 1 Patch]);
    
%     for ii = 1:Patch
%         ZFiller(:,:,ii) = lerp(SDF2+0.01*ii,SDF0,...
%                                 cos((1/Patch)*pi));
%     end
% end

% zero padding
SDF(:,:,1) = 0; SDF(:,:,end) = 0;
SDF(1,:,:) = 0; SDF(end,:,:) = 0;
SDF(:,1,:) = 0; SDF(:,end,:) = 0;

% Fem.Topology = SDF;

zq = linspace(0,Thickness,Layers);
[X,Y,Z] = meshgrid(xq,yq,zq);

% Fem.TopologyGridX = X; 
% Fem.TopologyGridY = Y; 
% Fem.TopologyGridZ = Z; 

end

function SDF = generateIsoCell(Fem,SDF0)
    vq = SDF0;
        
    if ~isempty(Fem.ReflectionPlane)
        RP = Fem.ReflectionPlane;
        if RP(1) == 1 && RP(2) == 1
            V = flip(vq,2); V = horzcat(V,vq);
            VV = flip(V,1); SDF = vertcat(VV,V);
        elseif RP(1) == 1 && RP(2) ~= 1
            V = flip(vq,2); V = horzcat(vq,V);
            SDF = V;
        elseif RP(1) == -1 && RP(2) ~= 1
            V = flip(vq,2); V = horzcat(V,vq);
            SDF = flip(V);
        elseif RP(1) ~= 1 && RP(2) == 1
            V = flip(vq,1); SDF = vertcat(V,vq);
        elseif RP(1) ~= 1 && RP(2) == -1
            V = flip(vq,1); SDF = vertcat(vq,V);
        elseif RP(1) == 0 && RP(2) == 0
            SDF = flip(vq);
        end
    else
        SDF = flip(vq);
    end
end