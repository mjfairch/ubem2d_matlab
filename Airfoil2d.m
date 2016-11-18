%
% Author:   Michael J. Fairchild (Princeton University)
% Date:     2016-09-22
%
classdef Airfoil2d < Body2d
    properties
        le;     % index of leading-edge node
        chord;  % chord length of airfoil
        xp;     % x-coordinate of pitch axis
        yp;     % y-coordinate of pitch axis
        puccw;  % true/false if pitch up is counterclockwise/clockwise
    end % properties
    
    methods
        % (x,y) = coordinates of airfoil nodes
        % le = index of leading-edge node
        % Points must be ordered counterclockwise starting at the trailing
        % edge and moving along the upper surface to the leading edge and
        % then back to the trailing edge along the lower surface.
        function obj = Airfoil2d(x,y,le,pp,puccw)
            obj@Body2d(x,y);    % call superclass constructor
            if (nargin > 0)
                if (nargin < 3)
                    error('Must specify x,y nodes and leading-edge index');
                end
                obj.le = le;
                chordVec = [x(1)-x(le),y(1)-y(le)];
                obj.chord = norm(chordVec);
                if (nargin < 4)
                    pp = 0;
                end
                obj.xp = x(le) + pp*chordVec(1);
                obj.yp = y(le) + pp*chordVec(2);
                if (nargin < 5)
                    if (x(le) > x(1))
                        obj.puccw = true;
                    elseif (x(le) < x(1))
                        obj.puccw = false;
                    else
                        error('Cannot determine pitch orientation');
                    end
                else
                    obj.puccw = puccw;
                end
            end
        end
        
%         function aoa = angleOfAttackLE(self,Ux,Uy,Vx,Vy,dalpdt,pp)
%             % TODO: Fix this method so that it works for all airfoil shapes
%             % and relative winds
%             error('Method not yet implemented');
%             % Inclination of the relative wind to the horizontal
%             [xp,yp] = self.chordline(pp);
%             [xle,yle] = self.leadingEdge();
%             U = [Ux,Uy];
%             V = [Vx + dalpdt*(yle-yp), Vy - dalpdt*(xle-xp)];
%             relWind = U - V;    % relative wind: flow as seen from airfoil
%             alpRelWind = atan2(relWind(2),relWind(1));
%             % Inclination of the chord line to the horizontal
%             [xte,yte] = self.trailingEdge();
%             chordDir = [xte-xle,yte-yle];
%             alpChord = atan2(chordDir(2),chordDir(1));
%             % Combine them
%             if (self.puccw)
%                 aoa = alpRelWind + alpChord;
%             else
%                 aoa = alpRelWind - alpChord;
%             end
%         end
        
        function [x,y] = getLeadingEdge(self)
            x = self.x(self.le);
            y = self.y(self.le);
        end
        
        function [x,y] = getTrailingEdge(self)
            x = self.x(end);
            y = self.y(end);
        end
        
        function theta = getTrailingEdgeBisector(self)
            v = .5*(self.panels{end}.tvec - self.panels{1}.tvec);
            theta = atan2(v(2),v(1));
        end
        
        function [x,y] = getPointOnChordLine(self,pp)
            x = self.x(self.le) + pp*(self.x(end) - self.x(self.le));
            y = self.y(self.le) + pp*(self.y(end) - self.y(self.le));
        end
        
        function [xp,yp] = getPitchAxis(self)
            xp = self.xp;
            yp = self.yp;
        end
        
        function setPitchAxis(self,xp,yp)
            self.xp = xp;
            self.yp = yp;
        end
        
        function [xp,yp] = setPitchAxisOnChordLine(self,pp)
            [xle,yle] = self.getLeadingEdge();
            [xte,yte] = self.getTrailingEdge();
            chordVec = [xte-xle,yte-yle];
            xp = self.x(self.le) + pp*chordVec(1);
            yp = self.y(self.le) + pp*chordVec(2);
            self.xp = xp;
            self.yp = yp;
        end
        
        function pitch(self,dtheta)
            self.glide(dtheta,0,0);
        end
        
        function heave(self,dy)
            self.glide(0,0,dy);
        end
        
        function surge(self,dx)
            self.glide(0,dx,0);
        end

        function glide(self,dtheta,dx,dy)
            if (nargin < 3)
                dx = 0;
                dy = 0;
            end
            if (self.puccw)
                glide@Body2d(self,dtheta,dx,dy,self.xp,self.yp);
                [xpnew,ypnew] = se2point(self.xp,self.yp,dx,dy,0);
            else
                glide@Body2d(self,-dtheta,dx,dy,self.xp,self.yp);
                [xpnew,ypnew] = se2point(self.xp,self.yp,dx,dy,0);
            end
            self.xp = xpnew;
            self.yp = ypnew;
        end
    end % methods
end