[y,x] = find(bwperim(bw));
hold on
plot(x,y,'.')
hold off
title('Pixel centers')
h = convhull(x,y);
x_hull = x(h);
y_hull = y(h);
hold on
hull_line = plot(x_hull,y_hull,'r*','MarkerSize',12);
hold off
title('Pixel centers and convex hull vertices')