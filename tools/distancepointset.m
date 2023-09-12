function d = distancepointset(PS1,PS2,R,varargin)
boundingBox      = [];
periodicBoundary = [];

assert(isa(1,'float') && R > 0, 'Radius must be a positive real scalar.');
assert(size(PS1,2) == size(PS1,2), 'PS1 and PS2 must be of same dimension.');

if isempty(boundingBox)
    B = boxhull([PS1;PS2]);
end

d = cell(size(PS1,1),1);
B1  = B(1); B2 = B(2); 
B12 = lerp(B1,B2,0.5);
B3  = B(3); B4 = B(4); 
B34 = lerp(B3,B4,0.5);

if size(PS1,2) == 3
    B5  = B(5); B6 = B(6);
    B56 = lerp(B5,B6,0.5);
end
    
for el = 1:size(PS1,1)   
    % if size(PS1, 2) == size(PS2, 2)
    %     dist = sqrt(sum((PS1(el,:) - PS2).^2, 2));
    % else
    %     error('Input matrices must have the same number of columns');
    % end
    if size(PS1,2) == 3
        dist = sqrt((PS1(el,1)-PS2(:,1)).^2 + (PS1(el,2)-PS2(:,2)).^2 +     ...
            (PS1(el,3)-PS2(:,3)).^2 );    
    else    
        dist = sqrt((PS1(el,1)-PS2(:,1)).^2 + (PS1(el,2)-PS2(:,2)).^2);
    end
    
    if ~isempty(periodicBoundary)
        Rp = periodicBoundary;
        dist2 = zeros(size(PS2,1),1);

        if Rp(1) == 1
            dist2(:,1) = sqrt((PS1(el,1)-(PS2(:,1) + B2-B1)).^2 + (PS1(el,2)-PS2(:,2)).^2);
        end

        if Rp(1) == 1/2
            PS = -(PS2(:,1)-B12) + B12;
            dist2(:,2) = sqrt((PS1(el,1)-PS).^2 + (PS1(el,2)-PS2(:,2)).^2);
        end

        if Rp(2) == 1/2
            PS = -(PS2(:,2)-B34) + B34;
            dist2(:,3) = sqrt((PS1(el,1)-PS2(:,1)).^2 + (PS1(el,2)-PS).^2);
        end

        dist = min([dist,dist2],[],2);
    end
    % if ~isempty(periodicBoundary)
    %     Rp = periodicBoundary;
    %     if Rp(1) == 1
    %          dist2 = sqrt((PS1(el,1)-(PS2(:,1) + B2-B1)).^2 +...
    %          (PS1(el,2)-PS2(:,2)).^2);
    %          dist = (min([dist,dist2].'))';
    %     end
    %     if Rp(1) == 1/2
    %         PS = -(PS2(:,1)-B12) + B12;
    %         dist2 = sqrt((PS1(el,1)-PS).^2 +...
    %          (PS1(el,2)-PS2(:,2)).^2);
    %          dist = (min([dist,dist2].'))';
    %     end
    %     if Rp(2) == 1/2
    %         PS = -(PS2(:,2)-B34) + B34;
    %         dist2 = sqrt((PS1(el,1)-PS2(:,1)).^2 +...
    %          (PS1(el,2)-PS).^2);
    %          dist = (min([dist,dist2].'))';
    %     end
    % end

    [I,J] = find(dist<=R);   
    d{el} = [I,J+(el-1),dist(I)];
end

% matrix of indices and distance value
d = cell2mat(d); 
end