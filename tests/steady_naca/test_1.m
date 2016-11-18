addpath('../../');
if (~exist('data','dir'))
    mkdir('data');
end

% Define test cases
codes = {'0012','2412','4412'}; % NACA 4-digit series
angles = [0,4,8];               % angles of attack (degrees)
pp = .25;                       % pitch axis/chord (0=LE, 1=TE)

% Store list of test cases for consumption by plotting script
fid = fopen('data/cases.txt','w');
for i=1:length(codes)
    fprintf(fid,'%04d ',str2num(codes{i}));
end
fprintf(fid,'\n');
for i=1:length(angles)
    fprintf(fid,'%g ',angles(i));
end
fprintf(fid,'\n');
fclose(fid);

% Run all test cases and store results
n = 50;
for i=1:length(codes)
    % Create body and save to disk
    code = codes{i};
    foil = foil_naca4(code,n,true);
    foil.setPitchAxisOnChordLine(pp);
    xmid = foil.getMidpoints();
    fid = fopen(sprintf('data/%s_body.txt',code),'w');
    fprintf(fid,'%13d %13d\n',foil.n,foil.le);
    for j=1:foil.n
        fprintf(fid,'%13.6g %13.6g\n',foil.x(j),foil.y(j));
    end
    fclose(fid);
    
    % Solve for flow past the body at various angles of attack
    for j = 1:length(angles)
        alp = angles(j);
        fprintf('%s, %g degrees\n',code,alp);
        % Solve steady-flow problem
        Uinf = [cos(alp*pi/180),sin(alp*pi/180)];
        [Cp,gam] = slvs2dhs(foil,Uinf);
        circ = gam*foil.perim;
        [xp,yp] = foil.getPitchAxis();
        [CFx,CFy,Cm] = forcemoment2d(foil,Cp,foil.chord,xp,yp,foil.puccw);
        [Cd,Cl] = aerocoef2d(CFx,CFy,Uinf);

        % Save pressure distribution to disk
        fid = fopen(sprintf('data/%s_%g_cp.txt',code,alp),'w');
        for k=1:foil.getNumberOfPanels()
            fprintf(fid,'%13.6g %13.6g\n',xmid(k),Cp(k));
        end
        fclose(fid);
        % Save aerodynamic coefficients to disk
        fid = fopen(sprintf('data/%s_%g_aero.txt',code,alp),'w');
        fprintf(fid,'%13.6g %13.6g %13.6g %13.6g %13.6g\n',circ,Cd,Cl,Cm,pp);
        fclose(fid);
    end
end