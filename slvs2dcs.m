% Solve steady two-dimensional flow problem by using the
% constant-strength source-panel method
% 
% Author: Michael J. Fairchild (Princeton University)
% Date: 2016-09-22
%
function [Cp,sigma] = slvs2dcs(body,An)
if (length(body) > 1)
    error('Multiple bodies not yet supported'); % TODO: implement it!
end

% Generate influence matrices if not provided
if (nargin < 2)
    [An,At] = inflmat2d(body.panels);
end
n = length(body.panels);

% Set up the right-hand side (for uniform rightward flow of unit speed)
Uinf = [1,0];
rhs = zeros(n,1);
for i=1:n
    rhs(i) = -dot(Uinf,body.panels{i}.nvec);
end

% Solve for panel source strengths
sigma = An\rhs;

% Compute the flow velocity at panel midpoints
q = At*sigma;
for i = 1:n
    q(i) = q(i) + dot(Uinf,body.panels{i}.tvec);
end

% Compute the pressure distribution with the steady Bernoulli equation
Cp = 1 - (q/norm(Uinf)).^2;