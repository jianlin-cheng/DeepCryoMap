function T = feretProperties(T)
% Copyright 2017-2018 The MathWorks, Inc.

maxfd = zeros(height(T),1);
maxfd_endpoints = cell(height(T),1);
maxfd_orientation = zeros(height(T),1);

minfd = zeros(height(T),1);
minfd_triangle_points = cell(height(T),1);
minfd_orientation = zeros(height(T),1);

minod = zeros(height(T),1);
minod_lower_points = cell(height(T),1);
minod_upper_points = cell(height(T),1);

minbb = cell(height(T),1);
minbb_a = zeros(height(T),1);

for k = 1:height(T)
    pixels = T.PixelList{k};
    V = pixelHull(pixels,'diamond');
    pairs = antipodalPairs(V);
    [maxfd(k),maxfd_endpoints{k}] = maxFeretDiameter(V,pairs);
    points = maxfd_endpoints{k};
    e = points(2,:) - points(1,:);
    maxfd_orientation(k) = atan2d(e(2),e(1));
    
    [minfd(k),minfd_triangle_points{k}] = minFeretDiameter(V,pairs);
    points = minfd_triangle_points{k};
    e = points(2,:) - points(1,:);
    thetad = atan2d(e(2),e(1));
    minfd_orientation(k) = mod(thetad + 180 + 90,360) - 180;
    
    [minod(k),minod_lower_points{k},minod_upper_points{k}] = ...
        feretDiameter(V,maxfd_orientation(k)+90);
    
    [minbb{k},minbb_a(k)] = minAreaBoundingBox(V,pairs);
end

T.MaxFeretDiameter = maxfd;
T.MaxFeretDiameterEndpoints = maxfd_endpoints;
T.MaxFeretDiameterOrientation = maxfd_orientation;
T.MinFeretDiameter = minfd;
T.MinFeretDiameterTrianglePoints = minfd_triangle_points;
T.MinFeretDiameterOrientation = minfd_orientation;
T.OrthogonalDiameter = minod;
T.OrthogonalDiameterLowerPoints = minod_lower_points;
T.OrthogonalDiameterUpperPoints = minod_upper_points;
T.MinAreaBoundingBox = minbb;
T.MinAreaBoundingBoxArea = minbb_a;
end

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

function h = triangleHeight(P1,P2,P3)
% Copyright 2017-2018 The MathWorks, Inc.

h = 2 * abs(signedTriangleArea(P1,P2,P3)) / norm(P1 - P2);
end

function area = signedTriangleArea(A,B,C)
% Copyright 2017-2018 The MathWorks, Inc.

area = ( (B(1) - A(1)) * (C(2) - A(2)) - ...
    (B(2) - A(2)) * (C(1) - A(1)) ) / 2;
end

function [d,V1,V2] = feretDiameter(V,theta)
% Copyright 2017-2018 The MathWorks, Inc.

% Rotate points so that the direction of interest is vertical.

alpha = 90 - theta;

ca = cosd(alpha);
sa = sind(alpha);
R = [ca -sa; sa ca];

% Vr = (R * V')';
Vr = V * R';

y = Vr(:,2);
ymin = min(y,[],1);
ymax = max(y,[],1);

d = ymax - ymin;

if nargout > 1
    V1 = V(y == ymin,:);
    V2 = V(y == ymax,:);
end
end