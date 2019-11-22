function vert = posToVert(posn,closed)
% Recover vertices from a rectangular position vector ([xmin ymin width height])
%
% Many functions in MATLAB and in Toolboxes specify rectangular regions in
% terms of "positions": [xmin ymin width height]. Many others require
% specification of the vertices of bounding regions. This function, and its
% companion |vertToPos|, facilitates easy conversion between them.
%
% SYNTAX:
% vert = posToVert(posn)
% vert = posToVert(posn,closed)
%
% INPUTS:
% posn: a 4-element vector indicating the position of a rectangular region,
%       as [xmin ymin width height]
%
% closed (optional): indicates whether the returned vertices should be
%       "closed", meaning that the last vertex will be repeated. Logical
%       true/false. (Default: false);
%
% OUTPUTS: 
% vert: a 4x2 (if 'open') or 5x2 (if 'closed') vector of the [x,y] vertices 
%       extracted from the region.
%
% EXAMPLES:
% %1) RECTANGLES and PATCHES
% % Rectangles are specified using "position," [x y w h]; patches are
% % specified using vertices. This example demonstrates that once either is
% % specified, you can calculate the other.
%
% t = 0:pi/32:6;
% plot(t,sin(t),'r-',t,cos(t),'b--')
% legend({'Sine','Cosine'});
% % FOR DEMONSTRATION, create two rectangles. The first we
% %   specify by [x y w h]:
% posn = [pi/4-0.25 sin(pi/4)-0.1 0.5 0.2];
% xy = posToVert(posn);
% rectangle('position',posn,...
% 	'linewidth',2,'edgecolor','r');
% h = patch(xy(:,1),xy(:,2),'c',...
% 	'facealpha',0.5,'linestyle','--','linewidth',2);
% % The second, we specify using vertices of polygon:
% xy = [5*pi/4-0.25 5*pi/4+0.25 5*pi/4+0.25 5*pi/4-0.25;
% 	sin(5*pi/4)-0.1 sin(5*pi/4)-0.1 sin(5*pi/4)+0.1 sin(5*pi/4)+0.1]';
% posn = vertToPos(xy);
% rectangle('Position',posn,...
% 	'linewidth',2,'edgecolor','r');
% patch(xy(:,1),xy(:,2),'y',...
% 	'facealpha',0.5,'linestyle','--','linewidth',2);
% title({'Rectangles are specified with position vectors';...
% 	'Patches are specified by vertices'})
%
% %2) IMRECT 
% % Note that this example uses Image Processing Toolbox functionality
% I = imread('eight.tif');
% subplot(1,2,1);
% imshow(I)
% t = title({'Drag an IMRECT around a quarter'; 'Double-click to continue.'});
% h = imrect; %Draw freehand rectangle
% position = wait(h);
% set(t,'string','Abracadabra...');
% vert = posToVert(position);
% J = regionfill(I,vert(:,1),vert(:,2));
% subplot(1,2,2);
% imshow(J);
% title('It''s magic!')
%
% Brett Shoelson, PhD.
% brett.shoelson@mathworks.com
% 12/08/2014
%
% See also: vertToPos, bbox2points

% Copyright 2014 The MathWorks, Inc.

narginchk(1,2)
if nargin == 1
	closed = false;
end
if numel(posn) ~= 4
	error('posToVert: posn must be a 4-element vector of the form [xmin ymin width height].')
end
validTypes = {'single','double','int16','int32', 'uint8', 'uint16', 'uint32'};
validateattributes(posn,validTypes,{'nonsparse','real','vector','nonnan'}, ...
    mfilename,'posn',1);

xs = [posn(1) posn(1)+posn(3) posn(1)+posn(3) posn(1)];
ys = [posn(2) posn(2) posn(2)+posn(4) posn(2)+posn(4)];
if closed
	xs = [xs xs(1)];
	ys = [ys ys(1)];
end
vert = [xs' ys'];
