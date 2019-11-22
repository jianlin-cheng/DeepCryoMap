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

function [d,end_points] = maxFeretDiameter(P,antipodal_pairs)
% Copyright 2017-2018 The MathWorks, Inc.

if nargin < 2
    antipodal_pairs = antipodalPairs(P);
end

v = P(antipodal_pairs(:,1),:) - P(antipodal_pairs(:,2),:);
D = hypot(v(:,1),v(:,2));
[d,idx] = max(D,[],1);
P1 = P(antipodal_pairs(idx,1),:);
P2 = P(antipodal_pairs(idx,2),:);

end_points = [P1 ; P2];
end

