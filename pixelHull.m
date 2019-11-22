function V = pixelHull(P,type)

if nargin < 2
    type = 24;
end

if islogical(P)
    P = bwperim(P);
    [i,j] = find(P);
    P = [j i];
end

if strcmp(type,'square')
    offsets = [ ...
         0.5  -0.5
         0.5   0.5
        -0.5   0.5
        -0.5  -0.5 ];

elseif strcmp(type,'diamond')
    offsets = [ ...
         0.5  0
         0    0.5
        -0.5  0
         0   -0.5 ];

else
    % type is number of angles for sampling a circle of diameter 1.
    thetad = linspace(0,360,type+1)';
    thetad(end) = [];

    offsets = 0.5*[cosd(thetad) sind(thetad)];
end

offsets = offsets';
offsets = reshape(offsets,1,2,[]);

Q = P + offsets;
R = permute(Q,[1 3 2]);
S = reshape(R,[],2);

k = convhull(S,'Simplify',true);
V = S(k,:);
end