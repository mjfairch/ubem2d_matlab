addpath('../../');
if (~exist('data','dir'))
    mkdir('data');
end

% Create body and panels
n = 50;
rad = 1;
t = linspace(0,2*pi,n+1);
x = rad*cos(t);
y = rad*sin(t);
cyl = Body2d(x,y);
n = cyl.getNumberOfPanels();

% Save body to disk
fid = fopen('data/cylinder_body.txt','w');
for i=1:length(t)
    fprintf(fid,'%13.6g %13.6g\n',x(i),y(i));
end
fclose(fid);

% Solve steady-flow problem
Cp = slvs2dcs(cyl);

% Save pressure distribution to disk
[xmid,ymid] = cyl.getMidpoints();
fid = fopen('data/cylinder_cp.txt','w');
for i=1:n
    fprintf(fid,'%13.6g %13.6g\n',atan2(ymid(i),xmid(i)),Cp(i));
end
fclose(fid);

fprintf('Done with cylinder test\n');