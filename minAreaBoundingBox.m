function [bb,A] = minAreaBoundingBox(V,antipodal_pairs)
% Copyright 2017-2018 The MathWorks, Inc.

if nargin < 2
    antipodal_pairs = antipodalPairs(V);
end

n = size(antipodal_pairs,1);
p = antipodal_pairs(:,1);
q = antipodal_pairs(:,2);

A = Inf;
thetad = [];

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
        % Points pt1, pt2, and pt3 are possibly on the minimum-area
        % bounding box, with points pt1 and pt2 forming an edge coincident with
        % the bounding box. Call the height of the triangle with base
        % pt1-pt2 the height of the bounding box, h.
        
        h = triangleHeight(pt1,pt2,pt3);
        
        pt1pt2_direction = atan2d(pt2(2)-pt1(2),pt2(1)-pt1(1));
        
        w = feretDiameter(V,pt1pt2_direction);
        
        A_k = h*w;
        if (A_k < A)
            A = A_k;
            thetad = pt1pt2_direction;
        end
    end
end

% Rotate all the points so that pt1-pt2 for the minimum bounding box points
% straight up.

r = 90 - thetad;
cr = cosd(r);
sr = sind(r);
R = [cr -sr; sr cr];

Vr = V * R';

xr = Vr(:,1);
yr = Vr(:,2);
xmin = min(xr);
xmax = max(xr);
ymin = min(yr);
ymax = max(yr);

bb = [ ...
    xmin ymin
    xmax ymin
    xmax ymax
    xmin ymax
    xmin ymin ];

% Rotate the bounding box points back.
bb = bb * R;
end

