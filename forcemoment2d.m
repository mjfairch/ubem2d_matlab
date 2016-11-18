%
% Author: Michael J. Fairchild (Princeton University)
% Date: 2016-09-25
%
% (x0,y0) = point about which to compute moments.  A counterclockwise
% moment is positive/negative according as to whether ccw is true/false.
function [CFx,CFy,Cm] = forcemoment2d(body,Cp,CharLen,x0,y0,ccw)
n = body.getNumberOfPanels();
Fx = zeros(1,n);
Fy = zeros(1,n);
for i=1:n
    Fx(i) = -Cp(i)*body.panels{i}.len*body.panels{i}.nvec(1);
    Fy(i) = -Cp(i)*body.panels{i}.len*body.panels{i}.nvec(2);
end
CFx = sum(Fx)/CharLen;
CFy = sum(Fy)/CharLen;
if (nargout > 2)
    [xmid,ymid] = body.getMidpoints();
    Cm = sum((Fy.*(xmid-x0) - Fx.*(ymid-y0)))/CharLen^2;    % CCW positive
    if (~ccw)
        Cm = -Cm;
    end
end