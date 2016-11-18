classdef UBEMStepper2d < handle
    properties
        % ========================================
        % Stepper parameters
        % ========================================
        % Reference point at "infinity" and computation of potential
        xref = -10; % x coordinate of point at infinity (as fraction of chord)
        yref = 0;   % y coordinate of point at infinity (as fraction of chord)
        nref = 20;  % number of panels from reference point to leading edge
        % Convergence parameters
        maxiters = 200; % max iterations for wake panel convergence
        tol = 1e-6;     % wake panel convergence tolerance
        quadtol = 1e-6; % truncation level for quadratic Kutta term
        maxerr = 1e-5;  % max error for Kutta or Neumann error norms
        % Wake modeling
        wakep_free = true;  % use free wake panel? (false: bisect trailing edge)
        wake_body = true;   % should body influence wake advection?
        wake_self = true;   % should wake influence its own advection?
        wake_rmin = 0;      % regularization: Rankine vortex below this distance

        % ========================================
        % State
        % ========================================
        foil;   % Airfoil2d object
        An;     % An(i,j) = norm. vel. at panel i midpt. due to src. panel j.
        At;     % At(i,j) = tang. vel. at panel i midpt. due to src. panel j.
        Bn;     % Bn(i,j) = norm. vel. at panel i midpt. due to vtx. panel j.
        Bt;     % Bt(i,j) = tang. vel. at panel i midpt. due to vtx. panel j.
        Uinf;   % 1x2 vector with x,y components of onset flow
        steps;  % the number of steps that have been taken
        gamma;  % circulation per unit length
        circb;  % bound circulation
        circt;  % total circulation (bound + wake)
        pot;    % potential at panel midpoints
        wpan;   % Panel2d object representing the wake panel
        wake;   % Wake2d object encoding the wake vortices
    end % properties
    
    methods
        % Constructor
        function obj = UBEMStepper2d(foil,Uinf)
            obj.steps = 0;
            obj.wake = Wake2d();
            if (nargin > 0)
                obj.foil = foil;
                obj.Uinf = Uinf;
            end
        end

        % Advance the simulation by one time step
        function [Cp,xp,yp] = step(self,dalp,dx,dy,dt)
            if (self.steps == 0)
                if (nargin < 2)
                    dalp = 0;
                    dx = 0;
                    dy = 0;
                end
                [Cp,phi,gam] = self.initialStep(dalp,dx,dy);
            else
                [Cp,phi,gam] = self.unsteadyStep(dalp,dx,dy,dt);
            end
            % Update bound and total circulation
            self.circb = gam*self.foil.perim;
            self.circt = self.circb + sum(self.wake.nu);
            % Prepare for the next step
            self.gamma = gam;
            self.pot = phi;
            self.steps = self.steps + 1;
            [xp,yp] = self.foil.getPitchAxis();
        end

        % Solve the steady-flow problem for the initial step
        function [Cp,phi,gam] = initialStep(self,dalp,dx,dy)
            % Move the body based on the kinematic arguments
            if (nargin > 1)
                self.foil.glide(dalp,dx,dy);
            end
            % Compute influence matrices
            [self.An,self.At,self.Bn,self.Bt] = inflmat2d(self.foil.panels);
            % Solve initial steady-flow problem
            [Cp,gam,sigma] = slvs2dhs(self.foil,self.Uinf,...
                self.An,self.At,self.Bn,self.Bt);
            % Compute potential
            qt = self.computeFlow(sigma,gam);
            phi = self.computePotential(qt,sigma,gam);
            % Initial guess for wake panel
            delk = self.foil.perim/self.foil.getNumberOfPanels();
            thk = self.foil.getTrailingEdgeBisector();
            self.updateWakePanel(delk,thk);
        end

        % Make a single unsteady step given kinematic changes
        function [Cp,phi,gamk] = unsteadyStep(self,dalp,dx,dy,dt)
            npan = self.foil.getNumberOfPanels();
            % Move the body based on the kinematic arguments
            self.foil.glide(dalp,dx,dy);
            % Solve implicit Kutta condition by iteration
            [xp,yp] = self.foil.getPitchAxis();
            Vn = self.surfaceNormalVelocity(xp,yp,dalp,dx,dy,dt);
            [gamk,sigk,delk,Uwk,Vwk] = self.solveImplicitKutta(Vn,dt);
            % Kelvin ciruclation theorem gives wake panel circulation
            L = self.foil.perim;
            gamwk = (L/delk)*(self.gamma-gamk);
            % Flow at panel midpoints
            [qt,qn] = self.computeFlow(sigk,gamk,gamwk);
            q = sqrt(qt.*qt + qn.*qn);
            % Check the boundary conditions and the Kutta condition
            errNeumann = norm(qn-Vn);
            errKutta = abs(q(1)^2-q(npan)^2-2*L*(gamk-self.gamma)/dt);
            if (errNeumann > self.maxerr)
                error('Neumann error: %g\n',errNeumann);
            end
            if (errKutta > self.maxerr)
                error('Kutta error: %g\n',errKutta);
            end
            % Potential at panel midpoints
            phi = self.computePotential(qt,sigk,gamk,gamwk);
            dphidt = (phi-self.pot)/dt;
            % Pressure distribution via unsteady Bernoulli equation
            spdinf = norm(self.Uinf);
            Cp = 1 - (q.*q + 2*dphidt)/spdinf^2;
            % Detach the wake panel and advect the wake
            self.detachWakePanel(Uwk,Vwk,gamwk,delk,dt);
            self.advectWake(sigk,gamk,dt);
        end

        % Iteratively solve the nonlinear, implicit Kutta condition
        function [gamk,sigk,delk,Uwk,Vwk] = solveImplicitKutta(self,Vn,dt)
            % Compute flow at panel midpoints due to onset flow and wake
            [Uinft,Uinfn] = self.flowOnset();
            [Wvt,Wvn] = self.flowWake();
            % Prepare for iteration
            L = self.foil.perim;
            npan = self.foil.getNumberOfPanels();
            delk = self.wpan.len;
            if (self.wakep_free)
                thk = self.wpan.theta;
            else
                thk = self.foil.trailingEdgeBisector();
            end
            converged = false;
            % Wake panel iteration loop
            for iter=1:self.maxiters
                % Update wake panel geometry
                self.updateWakePanel(delk,thk);
                % Flow components at panel midpoints due to wake panel
                [Wpt,Wpn] = self.flowWakePanel(1);
                % Solve for sigk in terms of unknown gamk
                bk = (L/delk)*Wpn'-sum(self.Bn,2);
                ck = -Uinfn' - (L/delk)*self.gamma*Wpn' - Wvn' + Vn';
                xxk = self.An\bk;
                yyk = self.An\ck;
                % Use the unsteady, quadratic Kutta condition to determine gamk
                alpha1 = dot(self.At(1,:),xxk) + sum(self.Bt(1,:)) - ...
                    (L/delk)*Wpt(1);
                beta1 = dot(self.At(1,:),yyk) + (L/delk)*self.gamma*Wpt(1) + ...
                    Wvt(1) + Uinft(1);
                alphaN = dot(self.At(npan,:),xxk) + sum(self.Bt(npan,:)) - ...
                    (L/delk)*Wpt(npan);
                betaN = dot(self.At(npan,:),yyk) + ...
                    (L/delk)*self.gamma*Wpt(npan) + Wvt(npan) + Uinft(npan);
                zeta = alpha1^2 - alphaN^2;
                eta = 2*(alpha1*beta1 - alphaN*betaN - L/dt);
                chi = beta1^2 - betaN^2 + 2*L*self.gamma/dt + ...
                    Vn(1)^2 - Vn(npan)^2;
                if (eta^2 - 4*zeta*chi < 0)
                    error('Negative discriminant');
                end
                if (abs(zeta) < self.quadtol)
                    % If leading coefficient of quadratic Kutta equation is 
                    % nearly zero, treat it instead as a linear equation.  
                    % This happens, for example, when the airfoil is 
                    % symmetric and the bisecting wake-panel model is used.
                    fprintf('Ignoring almost-zero quadratic Kutta term\n');
                    gamk = -chi/eta;
                else
                    % Quadratic formula
                    gamk = (-eta - sqrt(eta^2 - 4*zeta*chi))/(2*zeta);
                end
                sigk = gamk*xxk + yyk;
                % Compute resulting flow at wake panel midpoint
                Uwk = self.Uinf(1);
                Vwk = self.Uinf(2);
                for j=1:npan
                    [u,v,U,V] = pvi2d(self.foil.panels{j},...
                        self.wpan.xmid,self.wpan.ymid);
                    Uwk = Uwk + sigk(j)*u + gamk*U;
                    Vwk = Vwk + sigk(j)*v + gamk*V;
                end
                [u,v] = self.wake.influence(self.wpan.xmid,self.wpan.ymid);
                Uwk = Uwk + u;
                Vwk = Vwk + v;
                % Update wake panel length and angle
                delk = sqrt(Uwk^2+Vwk^2)*dt;
                if (self.wakep_free)
                    thk = atan2(Vwk,Uwk);
                end
                % Check for convergence
                if (iter > 1 && norm([Uwk-Uwk0,Vwk-Vwk0]) < self.tol)
                    converged = true;
                    break;
                end
                % Prepare for the next iteration
                Uwk0 = Uwk;
                Vwk0 = Vwk;
            end
            if (~converged)
                error('Failed to converge during step %d',self.steps);
            end
            self.updateWakePanel(delk,thk);
        end
        
        % Recompute wake panel given its length and orientation
        function updateWakePanel(self,delk,thk)
            [x1,y1] = self.foil.getTrailingEdge();
            x2 = x1 + delk*cos(thk);
            y2 = y1 + delk*sin(thk);
            self.wpan = Panel2d(x1,y1,x2,y2); % ccw/cw orientation irrelevant
        end
        
        % Shed the wake panel as a point vortex
        function detachWakePanel(self,Uwk,Vwk,gamwk,delk,dt)
            % Detach wake panel (see equation (14) in Basu & Hancock)
            shedx = self.wpan.xmid + Uwk*dt;
            shedy = self.wpan.ymid + Vwk*dt;
            self.wake.addVortex(gamwk*delk, shedx, shedy);
        end
        
        % Determine net flow at each wake vortex and advect the wake
        function advectWake(self,sigk,gamk,dt)
            % Compute net flow velocity at all wake vortices
            npan = self.foil.getNumberOfPanels();
            nvort = self.wake.size();
            vx = self.Uinf(1)*ones(1,nvort);
            vy = self.Uinf(2)*ones(1,nvort);
            if (self.wake_body)
                for i = 1:npan
                    [u,v,U,V] = pvi2d(self.foil.panels{i},self.wake.x,...
                        self.wake.y);
                    vx = vx + sigk(i)*u + gamk*U;
                    vy = vy + sigk(i)*v + gamk*V;
                end
            end
            if (self.wake_self)
                [u,v] = self.wake.selfInfluence();
                vx = vx + u;
                vy = vy + v;
            end
            % Now advect the wake by the given time step
            self.wake.advect(vx,vy,dt);
        end
        
        % Compute velocity of panel midpoints due to kinematics
        function Vn = surfaceNormalVelocity(self,xp,yp,dalp,dx,dy,dt)
            om = dalp/dt;
            vx = dx/dt;
            vy = dy/dt;
            [xmid,ymid] = self.foil.getMidpoints();
            Vs = [vx + om*(ymid-yp); vy - om*(xmid-xp)];
            npan = self.foil.getNumberOfPanels();
            Vn = zeros(1,npan);
            for i=1:npan
                Vn(i) = dot(Vs(:,i),self.foil.panels{i}.nvec);
            end
        end
        
        % Compute net flow at panel midpoints after unsteady step
        function [qt,qn] = computeFlow(self,sigk,gamk,gamwk)
            [Uinft,Uinfn] = flowOnset(self);
            [Pant,Pann] = flowPanels(self,sigk,gamk);
            qt = Pant + Uinft;
            qn = Pann + Uinfn;
            if (self.steps > 0)
                [Wpt,Wpn] = self.flowWakePanel(gamwk);
                [Wvt,Wvn] = self.flowWake();
                qt = qt + Wpt + Wvt;
                qn = qn + Wpn + Wvn;
            end
        end
        
        % Compute flow at panel midpoints due to onset flow
        function [Uinft,Uinfn] = flowOnset(self)
            npan = self.foil.getNumberOfPanels();
            Uinft = zeros(1,npan);
            Uinfn = zeros(1,npan);
            for i=1:npan
                Uinft(i) = dot(self.Uinf,self.foil.panels{i}.tvec);
                Uinfn(i) = dot(self.Uinf,self.foil.panels{i}.nvec);
            end
        end
        
        % Compute flow at panel midpoints due to body panels
        function [Pant,Pann] = flowPanels(self,sigk,gamk)
            Pant = (self.At*sigk + gamk*sum(self.Bt,2))';
            Pann = (self.An*sigk + gamk*sum(self.Bn,2))';
        end
        
        % Compute flow at panel midpoints due to wake panel
        function [Wpt,Wpn] = flowWakePanel(self,gamwp)
            [xmid,ymid] = self.foil.getMidpoints();
            [~,~,U,V] = pvi2d(self.wpan,xmid,ymid);
            npan = self.foil.getNumberOfPanels();
            Wpt = zeros(1,npan);
            Wpn = zeros(1,npan);
            for i=1:npan
                Wpt(i) = gamwp*dot([U(i),V(i)],self.foil.panels{i}.tvec);
                Wpn(i) = gamwp*dot([U(i),V(i)],self.foil.panels{i}.nvec);
            end
        end
        
        % Compute flow at panel midpoints due to wake vortices
        function [Wvt,Wvn] = flowWake(self)
            [xmid,ymid] = self.foil.getMidpoints();
            [vx,vy] = self.wake.influence(xmid,ymid);
            npan = self.foil.getNumberOfPanels();
            Wvt = zeros(1,npan);
            Wvn = zeros(1,npan);
            for i=1:npan
                Wvt(i) = dot([vx(i),vy(i)],self.foil.panels{i}.tvec);
                Wvn(i) = dot([vx(i),vy(i)],self.foil.panels{i}.nvec);
            end
        end

        % Compute potential at panel midpoints by integrating the flow from
        % a far-upstream reference point to the panel midpoints
        function phi = computePotential(self,qt,sigk,gamk,gamwk)
            % Integrate the flow from the reference point to the leading edge
            [xle,yle] = self.foil.getLeadingEdge();
            xpp = linspace(self.xref,xle,self.nref+1);
            ypp = linspace(self.yref,yle,self.nref+1);
            uu = self.Uinf(1)*ones(1,self.nref);
            vv = self.Uinf(2)*ones(1,self.nref);
            npan = self.foil.getNumberOfPanels();
            for i = 1:npan
                [u,v,U,V] = pvi2d(self.foil.panels{i},xpp(1:end-1),...
                    ypp(1:end-1));
                uu = uu + sigk(i)*u + gamk*U;
                vv = vv + sigk(i)*v + gamk*V;
            end
            if (self.steps > 0)
                [U,V] = self.wake.influence(xpp(1:end-1),ypp(1:end-1));
                uu = uu + U;
                vv = vv + V;
                [~,~,U,V] = pvi2d(self.wpan,xpp(1:end-1),ypp(1:end-1));
                uu = uu + gamwk*U;
                vv = vv + gamwk*V;
            end
            le = self.foil.le;
            phi = zeros(1,npan+1);
            phi(le) = sum(uu.*diff(xpp) + vv.*diff(ypp));
            % Potential at upper surface nodes
            for i=le-1:-1:1
                phi(i) = phi(i+1) - qt(i)*self.foil.panels{i}.len;
            end
            % Potential at lower surface nodes
            for i=le+1:npan+1
                phi(i) = phi(i-1) + qt(i-1)*self.foil.panels{i-1}.len;
            end
            % Potential at midpoints: average over adjacent nodes
            phi = .5*(phi(1:end-1) + phi(2:end));
        end
    end % methods
end