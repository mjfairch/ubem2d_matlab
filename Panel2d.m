%
% Author:   Michael J. Fairchild (Princeton University)
% Date:     2016-09-22
classdef Panel2d < handle
    properties
        x1;     % (x1,y1) = coordinates of first panel corner
        y1;
        x2;     % (x2,y2) = coordinates of second panel corner
        y2;
        ccw;    % true/false if nvec rotated CCW/CW from tvec
        xmid;   % (xmid,ymid) = coordinates of panel midpoint
        ymid;
        len;    % length of panel
        tvec;   % 1x2 vector with x,y components of unit tangent vector
        nvec;   % 1x2 vector with x,y components of unit normal vector
        theta;  % angle between unit tangent vector and positive x axis
        beta;   % angle between unit normal vector and positive x axis
    end % properties
    
    methods
        % (x1,y1):  coordinates of first panel corner
        % (x2,y2):  coordinates of second panel corner
        % ccw:      true/false if normal vector should be rotated CCW/CW
        %           from tangent vector.  Default = false.
        function obj = Panel2d(x1,y1,x2,y2,ccw)
            if (nargin > 0)
                if (nargin < 5)
                    ccw = false;
                end
                obj.x1 = x1;
                obj.y1 = y1;
                obj.x2 = x2;
                obj.y2 = y2;
                obj.ccw = ccw;
                obj.update();
            end
        end
        
        function setFirstCorner(x1,y1)
            obj.x1 = x1;
            obj.y1 = y1;
            obj.update();
        end
        
        function setSecondCorner(x2,y2)
            obj.x2 = x2;
            obj.y2 = y2;
            obj.update();
        end
        
        function setOrientation(self,ccw)
            self.ccw = ccw;
            self.update();
        end
        
        function update(self)
            self.xmid = .5*(self.x1+self.x2);
            self.ymid = .5*(self.y1+self.y2);
            self.len = norm([self.x2-self.x1,self.y2-self.y1]);
            self.tvec = [self.x2-self.x1,self.y2-self.y1]/self.len;
            self.nvec = [self.tvec(2), -self.tvec(1)]; % clockwise rotation
            if (self.ccw)
                self.nvec = -self.nvec;
            end
            self.theta = atan2(self.tvec(2),self.tvec(1));
            self.beta = atan2(self.nvec(2),self.nvec(1));
        end
    end % methods
end