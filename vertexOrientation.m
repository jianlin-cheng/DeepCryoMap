function s = vertexOrientation(P0,P1,P2)
% vertexOrientation  Orientation of a vertex with respect to line segment.
%
%   s = vertexOrientation(P0,P1,P2) returns a positive number if P2 is to
%   the left of the line through P0 to P1. It returns 0 if P2 is on the
%   line. It returns a negative number if P2 is to the right of the line.
%
%   Stating it another way, a positive output corresponds to a
%   counterclockwise traversal from P0 to P1 to P2.
%
%   P0, P1, and P2 are two-element vectors containing (x,y) coordinates.
%
%   Reference: http://geomalgorithms.com/a01-_area.html, function isLeft()

% Steve Eddins


s = (P1(1) - P0(1)) * (P2(2) - P0(2)) - ...
    (P2(1) - P0(1)) * (P1(2) - P0(2));
end