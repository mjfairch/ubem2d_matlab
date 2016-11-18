% Panel Velocity Influence: returns the velocity induced by unit-strength
% source and vortex panels at a collection of points
%
% Input:    panel = panel object
%           x = vector of x coordinates at which to evaluate induced velocity
%           y = vector of y coordinates at which to evaluate induced velocity
%
% Output:   [u,v] = velocity due to unit source distribution along panel
%           [U,V] = velocity due to unit vortex distribution along panel
%
% Note:     (U,V) is related to (u,v) by: U = -v, V = u
%
% Example:  pan = Panel2d(0,0,1,1);
%           [u,v,U,V] = pvi2d(pan,[2,2],[2,-2]);
%
% Author:   Michael J. Fairchild (Princeton University)
% Date:     2016-09-22
function [u,v,U,V] = pvi2d(panel,x,y)
L = panel.len;
cth = cos(panel.theta);
sth = sin(panel.theta);
N = length(x);
if (length(y) ~= N)
    error('Size mismatch');
end
u = zeros(1,N);
v = zeros(1,N);
for i=1:N
    dx = x(i)-panel.x1;
    dy = y(i)-panel.y1;
    b = -2*(dx*cth + dy*sth);
    c = dx*dx + dy*dy;
    u(i) = u(i) + 1/(2*pi)*panelint2d(-cth,dx,b,c,L);
    v(i) = v(i) + 1/(2*pi)*panelint2d(-sth,dy,b,c,L);
end
U = -v;
V = u;