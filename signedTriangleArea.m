function area = signedTriangleArea(A,B,C)
% Copyright 2017-2018 The MathWorks, Inc.

area = ( (B(1) - A(1)) * (C(2) - A(2)) - ...
    (B(2) - A(2)) * (C(1) - A(1)) ) / 2;
end

