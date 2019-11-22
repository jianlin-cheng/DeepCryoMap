function h = triangleHeight(P1,P2,P3)
% Copyright 2017-2018 The MathWorks, Inc.

h = 2 * abs(signedTriangleArea(P1,P2,P3)) / norm(P1 - P2);
end