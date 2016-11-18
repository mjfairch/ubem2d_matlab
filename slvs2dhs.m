% Solve steady-state two-dimensional flow problem by using the Hess-Smith
% method, which employs constant-strength source and vortex panels.
%
% Input:
%   body : Body2d object
%   Uinf : 1x2 vector giving x,y components of onset flow
%
% Output:
%   Cp : vector of pressure coefficients, one entry for each panel
%   gamma : ciruclation (positive counterclockwise) per unit length
%   sigma : vector of source strengths, one entry for each panel
% 
% Author: Michael J. Fairchild (Princeton University)
% Date: 2016-09-22
%
function [Cp,gamma,sigma] = slvs2dhs(body,Uinf,An,At,Bn,Bt)
if (length(body) > 1)
    error('Multiple bodies not yet supported'); % TODO: implement it!
end
if (nargin < 2)
    Uinf = [1,0];
end
if (nargin > 2 && nargin < 6)
    error('Must specify all influence matrices');
end

% Generate influence matrices if not provided
if (nargin < 6)
    [An,At,Bn,Bt] = inflmat2d(body.panels);
end
n = length(body.panels);

% Set up the right-hand side
rhs = zeros(n+1,1);
for i=1:n
    rhs(i) = -dot(Uinf,body.panels{i}.nvec);
end
rhs(n+1) = -dot(Uinf,body.panels{1}.tvec + body.panels{n}.tvec);

% Set up the Hess-Smith matrix
A = zeros(n+1,n+1);
A(1:n,1:n) = An;
A(1:n,n+1) = sum(Bn,2);
A(n+1,1:n) = At(1,:) + At(n,:);
A(n+1,n+1) = sum(Bt(1,:)) + sum(Bt(n,:));

% Solve for source strengths and vortex strength
soln = A\rhs;
sigma = soln(1:n);
gamma = soln(n+1);

% Compute the flow velocity at panel midpoints
qt = At*sigma + gamma*sum(Bt,2);
qn = An*sigma + gamma*sum(Bn,2);
for i=1:n
    qt(i) = qt(i) + dot(Uinf,body.panels{i}.tvec);
    qn(i) = qn(i) + dot(Uinf,body.panels{i}.nvec);
end
q = sqrt(qt.*qt + qn.*qn);

% Compute the pressure distribution with the steady Bernoulli equation
Cp = 1 - (q/norm(Uinf)).^2;