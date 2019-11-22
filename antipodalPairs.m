
function pq = antipodalPairs(S)
% antipodalPairs Antipodal vertex pairs of simple, convex polygon.
%
%   pq = antipodalPairs(S) computes the antipodal vertex pairs of a simple,
%   convex polygon. S is a Px2 matrix of (x,y) vertex coordinates for the
%   polygon. S must be simple and convex without repeated vertices. It is
%   not checked for satisfying these conditions. S can either be closed or
%   not. The output, pq, is an Mx2 matrix representing pairs of vertices in
%   S. The coordinates of the k-th antipodal pair are S(pq(k,1),:) and
%   S(pq(k,2),:).
%
%   TERMINOLOGY
%
%   For a convex polygon, an antipodal pair of vertices is one where you
%   can draw distinct lines of support through each vertex such that the
%   lines of support are parallel.
%
%   A line of support is a line that goes through a polygon vertex such
%   that the interior of the polygon lies entirely on one side of the line.
%
%   EXAMPLE
%
%     Compute antipodal vertices of a polygon and plot the corresponding
%     line segments.
%
%       x = [0 0 1 3 5 4 0];
%       y = [0 1 4 5 4 1 0];
%       S = [x' y'];
%       pq = antipodalPairs(S);
%
%       plot(S(:,1),S(:,2))
%       hold on
%       for k = 1:size(pq,1)
%           xk = [S(pq(k,1),1) S(pq(k,2),1)];
%           yk = [S(pq(k,1),2) S(pq(k,2),2)];
%           plot(xk,yk,'LineStyle','--','Marker','o','Color',[0.7 0.7 0.7])
%       end
%       hold off
%       axis equal
%
%   ALGORITHM NOTES
%
%   This function uses the "ANTIPODAL PAIRS" algorithm, Preparata and
%   Shamos, Computational Geometry: An Introduction, Springer-Verlag, 1985,
%   p. 174.

%   Steve Eddins


n = size(S,1);

if isequal(S(1,:),S(n,:))
    % The input polygon is closed. Remove the duplicate vertex from the
    % end.
    S(n,:) = [];
    n = n - 1;
end

% The algorithm assumes the input vertices are in counterclockwise order.
% If the vertices are in clockwise order, reverse the vertices.
clockwise = simplePolygonOrientation(S) < 0;
if clockwise
    S = flipud(S);
end

% The following variables, including the two anonymous functions, are set
% up to follow the notation in the pseudocode on page 174 of Preparata and
% Shamos. p and q are indices (1-based) that identify vertices of S. p0 and
% q0 identify starting vertices for the algorithm. area(i,j,k) is the area
% of the triangle with the corresponding vertices from S: S(i,:), S(j,:),
% and S(k,:). next(p) returns the index of the next vertex of S.
%
% The initialization of p0 is missing from the Preparata and Shamos text.
area = @(i,j,k) signedTriangleArea(S(i,:),S(j,:),S(k,:));
next = @(i) mod(i,n) + 1; % mod((i-1) + 1,n) + 1
p = n;
p0 = next(p);
q = next(p);

% The list of antipodal vertices will be built up in the vectors pp and qq.
pp = zeros(0,1);
qq = zeros(0,1);

% ANTIPODAL PAIRS step 3.
while (area(p,next(p),next(q)) > area(p,next(p),q))
    q = next(q);
end
q0 = q;    % Step 4.

while (q ~= p0)    % Step 5.
    p = next(p);   % Step 6.
    % Step 7. (p,q) is an antipodal pair.
    pp = [pp ; p];
    qq = [qq ; q];

    % Step 8.
    while (area(p,next(p),next(q)) > area(p,next(p),q))
        q = next(q);    % Step 9.
        if ~isequal([p q],[q0,p0])
            % Step 10.
            pp = [pp ; p];
            qq = [qq ; q];
        else
            % This loop break is omitted from the Preparata and Shamos
            % text.
            break
        end
    end

    % Step 11. Check for parallel edges.
    if (area(p,next(p),next(q)) == area(p,next(p),q))
        if ~isequal([p q],[q0 n])
            % Step 12. (p,next(q)) is an antipodal pair.
            pp = [pp ; p];
            qq = [qq ; next(q)];
        else
            % This loop break is omitted from the Preparata and Shamos
            % text.
            break
        end
    end
end

if clockwise
    % Compensate for the flipping of the polygon vertices.
    pp = n + 1 - pp;
    qq = n + 1 - qq;
end

pq = [pp qq];
end