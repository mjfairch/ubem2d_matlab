%
% Author: Michael J. Fairchild (Princeton University)
% Date: 2016-09-27
%
function [xx,yy] = se2point(x,y,a,b,theta,x0,y0)
if (nargin < 6)
    x0 = 0;
    y0 = 0;
end
if (nargin < 5)
    theta = 0;
end
if (nargin > 2)
    cth = cos(theta);
    sth = sin(theta);
    xx = (x-x0)*cth - (y-y0)*sth + x0 + a;
    yy = (x-x0)*sth + (y-y0)*cth + y0 + b;
else
    xx = x;
    yy = y;
end