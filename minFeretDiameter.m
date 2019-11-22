function [d,triangle_points] = minFeretDiameter(V,antipodal_pairs)
% Copyright 2017-2018 The MathWorks, Inc.

if nargin < 2
    antipodal_pairs = antipodalPairs(V);
end

n = size(antipodal_pairs,1);
p = antipodal_pairs(:,1);
q = antipodal_pairs(:,2);

d = Inf;
triangle_points = [];

for k = 1:n
    if k == n
        k1 = 1;
    else
        k1 = k+1;
    end
    
    pt1 = [];
    pt2 = [];
    pt3 = [];
    
    if (p(k) ~= p(k1)) && (q(k) == q(k1))
        pt1 = V(p(k),:);
        pt2 = V(p(k1),:);
        pt3 = V(q(k),:);
        
    elseif (p(k) == p(k1)) && (q(k) ~= q(k1))
        pt1 = V(q(k),:);
        pt2 = V(q(k1),:);
        pt3 = V(p(k),:);
    end
    
    if ~isempty(pt1)
        % Points pt1, pt2, and pt3 form a possible minimum Feret diameter.
        % Points pt1 and pt2 form an edge parallel to caliper direction.
        % The Feret diameter orthogonal to the pt1-pt2 edge is the height
        % of the triangle with base pt1-pt2.
        
        d_k = triangleHeight(pt1,pt2,pt3);
        if d_k < d
            d = d_k;
            triangle_points = [pt1; pt2; pt3];
        end
    end
end
end    