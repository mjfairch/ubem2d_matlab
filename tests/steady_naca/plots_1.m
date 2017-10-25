if (~exist('plots','dir'))
    mkdir('plots');
end
addpath('../../') % for saveTightFigure
% Identify test cases
data = load('data/cases.txt');
nbodies = length(data(1,:));
nangles = length(data(2,:));
bodies{nbodies} = [];
for i=1:nbodies
    bodies{i} = sprintf('%04d',data(1,i));
end
angles = data(2,:);

for i=1:nbodies
    body = bodies{i};
    % Plot the airfoil section
    data = load(sprintf('data/%s_body.txt',body));
    figure;
    hold on;
    f = Airfoil2d(data(2:end,1),data(2:end,2),data(1,2));
    fill(f.x/f.chord,f.y/f.chord,[.9 .9 .9]);
    plot(f.x/f.chord,f.y/f.chord,'-ok','LineWidth',1.5);
    axis equal;
    grid on;
    box on;
    xlabel('$x/c$','Interpreter','latex');
    ylabel('$y/c$','Interpreter','latex');
    title(sprintf('NACA %s, %d panels',body,f.getNumberOfPanels()),...
        'Interpreter','latex');
    set(gca,'FontSize',14);
    saveTightFigure(sprintf('plots/%s_body.pdf',body));
    close;
    
    for j=1:nangles
        angle = angles(j);
        % Load pressure data
        data = load(sprintf('data/%s_%g_cp.txt',body,angle));
        xmid = f.getMidpoints();
        Cp = data(:,2);
        % Load aerodynamic data
        data = load(sprintf('data/%s_%g_aero.txt',body,angle));
        Cd = data(1,2);
        Cl = data(1,3);
        Cm = data(1,4);
        pp = data(1,5);
        figure;
        hold on;
        lep = f.le - 1;
        plot(xmid(1:lep)/f.chord,Cp(1:lep),'-ok','LineWidth',1.5,...
            'MarkerSize',12);
        plot(xmid(lep:end)/f.chord,Cp(lep:end),'-dr','LineWidth',1.5,...
            'MarkerSize',12);
        legend({'Upper surface','Lower surface'},'Interpreter','latex');
        set(gca,'Ydir','reverse');
        grid on;
        box on;
        xlabel('$x/c$','Interpreter','Latex');
        ylabel('$C_p$','Interpreter','latex');
        title([
            sprintf('NACA %s at %g degrees, %d panels\n',body,angle,f.n-1),...
            sprintf('$C_D$=%.3g, $C_L$=%.3g, $C_M$=%.3g (at x/c=%.3g)',...
            Cd,Cl,Cm,pp)],'Interpreter','latex');
        set(gca,'FontSize',14);
        saveTightFigure(sprintf('plots/%s_%g_cp.pdf',body,angle));
        close;
    end
end