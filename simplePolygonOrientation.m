function s = simplePolygonOrientation(V)
% simplePolygonOrientation  Determine vertex order for simple polygon.
%
%   s = simplePolygonOrientation(V) returns a positive number if the simple
%   polygon V is counterclockwise. It returns a negative number of the
%   polygon is clockwise. It returns 0 for degenerate cases. V is a Px2
%   matrix of (x,y) vertex coordinates.
%
%   Reference: http://geomalgorithms.com/a01-_area.html, function
%   orientation2D_Polygon()

% Steve Eddins


n = size(V,1);

if n < 3
    s = 0;
    return
end

% Find rightmost lowest vertext of the polygon.

x = V(:,1);
y = V(:,2);
ymin = min(y,[],1);
y_idx = find(y == ymin);
if isscalar(y_idx)
    idx = y_idx;
else
    [~,x_idx] = max(x(y_idx),[],1);
    idx = y_idx(x_idx(1));
end

% The polygon is counterclockwise if the edge leaving V(idx,:) is left of
% the entering edge.

if idx == 1
    s = vertexOrientation(V(n,:), V(1,:), V(2,:));
elseif idx == n
    s = vertexOrientation(V(n-1,:), V(n,:), V(1,:));
else
    s = vertexOrientation(V(idx-1,:), V(idx,:), V(idx+1,:));
end
end
