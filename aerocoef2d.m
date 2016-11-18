%
% Author:   Michael J. Fairchild (Princeton University)
% Date:     2016-09-22
%
% DragDir = 1x2 vector indicating the drag direction
% CCW: true/false if lift direction is rotated counterclockwise/clockwise from
% drag direction. (Default = true)
function [CD,CL] = aerocoef2d(CFx,CFy,DragDir,ccw)
if (nargin < 4)
    ccw = true;
end
if (nargin < 3)
    DragDir = [1,0];
end
% Normalize drag direction and rotate it to get the lift direction
DragDir = DragDir/norm(DragDir);
if (ccw)
    LiftDir = [-DragDir(2),DragDir(1)];
else
    LiftDir = [DragDir(2),-DragDir(1)];
end
CD = dot([CFx,CFy],DragDir);
CL = dot([CFx,CFy],LiftDir);