classdef Wake2d < handle
    properties
        nu;     % srength of vortices
        x;      % horizontal position of vortices
        y;      % vertical position of vortices
    end % properties
    
    methods
        function obj = Wake2d(nu,x,y)
            if (nargin > 1)
                n = length(nu);
                if (length(x) ~= n || length(x) ~= n)
                    error('Size mismatch');
                end
                obj.nu = nu;
                obj.x = x;
                obj.y = y;
            else
                obj.nu = [];
                obj.x = [];
                obj.y = [];
            end
        end
        
        function n = size(self)
            n = length(self.nu);
        end
        
        function addVortex(self,str,x0,y0)
            self.nu = [self.nu str];
            self.x = [self.x x0];
            self.y = [self.y y0];
        end
        
        function [str,x,y] = getVortices(self,minStr,maxStr)
            I = intersect(find(self.nu > minStr),find(self.nu < maxStr));
            str = self.nu(I);
            x = self.x(I);
            y = self.y(I);
        end
        
        function [str,x,y] = getPositiveVortices(self)
            [str,x,y] = self.getVortices(0,Inf);
        end
        
        function [str,x,y] = getNegativeVortices(self)
            [str,x,y] = self.getVortices(-Inf,0);
        end
        
        function [vx,vy] = influence(self,x,y)
            if (length(x) ~= length(y))
                error('Size mismatch');
            end
            n = length(x);
            vx = zeros(1,n);
            vy = zeros(1,n);
            for i=1:length(self.nu)
                dx = x - self.x(i);
                dy = y - self.y(i);
                dr2 = dx.*dx + dy.*dy;
                vx = vx - self.nu(i).*dy./(2*pi*dr2);
                vy = vy + self.nu(i).*dx./(2*pi*dr2);
            end
        end

        function [u,v] = selfInfluence(self)
            n = self.size();
            u = zeros(1,n);
            v = zeros(1,n);
            for i=1:n
                xi = self.x(i);
                yi = self.y(i);
                ui = 0;
                vi = 0;
                for j=1:n
                    if (j ~= i)
                        dx = xi - self.x(j);
                        dy = yi - self.y(j);
                        dr2 = dx*dx + dy*dy;
                        ui = ui - self.nu(j)*dy/(2*pi*dr2);
                        vi = vi + self.nu(j)*dx/(2*pi*dr2);
                    end
                end
                u(i) = ui;
                v(i) = vi;
            end
        end
        
        function advect(self,vx,vy,dt)
            self.x = self.x + vx*dt;
            self.y = self.y + vy*dt;
        end
        
        % vortexCores returns the vortex core strengths and centers.
        % A vortex core is delineated by a sequence of vortices whose strenghts 
        % are of constant sign (i.e. the boundary between two vortex cores is 
        % the place where the sequence of vortex strenghts changes sign).
        % The strength of the core is the sum of its constituent vortex 
        % strengths, and its center is defined by analogy with the 
        % center-of-mass formula, with vortex strength playing the role of 
        % mass.  Vortices of zero strength are omitted from consideration.
        %
        % Output:
        %   str     Vector of vortex core strengths
        %   cx      Vector of x coordinates of vortex cores
        %   cy      Vector of y coordinates of vortex cores
        %
        % Author: Michael J. Fairchild
        % Date: July, 2016
        %
        function [str,cx,cy] = vortexCores(self)
            % Eliminate "vortices" of zero strength
            in0 = find(self.nu ~= 0);
            mu = self.nu(in0);
            xx = self.x(in0);
            yy = self.y(in0);

            % Find indices where sign of vorticity changes
            s = diff(sign(mu));
            isc = find(s ~= 0);
            if (isempty(isc))
                str = sum(mu);
                cx = sum(mu.*xx)/str;
                cy = sum(mu.*yy)/str;
            else
                n = length(isc) + 1;    % # of vortex cores = # of sign changes + 1

                % Vortex core strengths and locations are defined analogously to to center of
                % mass, where core strength is the sum of the vortex strenghts in the core,
                % and the core locations are the weighted average of vortex positions.
                cx = zeros(1,n);
                cy = zeros(1,n);
                str = zeros(1,n);
                for i=1:n
                    if (i == 1)
                        i1 = 1;
                    else
                        i1 = isc(i-1)+1;
                    end
                    if (i == 1)
                        i2 = isc(i);
                    elseif (i == n)
                        i2 = length(mu);
                    else
                        i2 = isc(i);
                    end
                    str(i) = sum(mu(i1:i2));
                    cx(i) = sum(mu(i1:i2).*xx(i1:i2))/str(i);
                    cy(i) = sum(mu(i1:i2).*yy(i1:i2))/str(i);
                end
            end
        end
    end % methods
end