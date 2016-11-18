% 
% Author: Michael J. Fairchild (Princeton University)
% Date: 2016-09-22
%
function [An,At,Bn,Bt] = inflmat2d(panels)
% Compute panel-to-panel influence matrix
n = length(panels);
An = zeros(n,n);    % An(i,j) = norm. vel. at panel i midpt. due to src. panel j
At = zeros(n,n);    % At(i,j) = tang. vel. at panel i midpt. due to src. panel j
Bn = zeros(n,n);    % Bn(i,j) = norm. vel. at panel i midpt. due to vtx. panel j
Bt = zeros(n,n);    % Bn(i,j) = tang. vel. at panel i midpt. due to vtx. panel j
for i=1:n
    for j=1:n
        if (j == i)
            An(i,j) = .5;
            Bt(i,j) = .5;
        else
            [u,v,U,V] = pvi2d(panels{j},panels{i}.xmid,panels{i}.ymid);
            An(i,j) = dot([u,v],panels{i}.nvec);
            At(i,j) = dot([u,v],panels{i}.tvec);
            Bn(i,j) = dot([U,V],panels{i}.nvec);
            Bt(i,j) = dot([U,V],panels{i}.tvec);
        end
    end
end