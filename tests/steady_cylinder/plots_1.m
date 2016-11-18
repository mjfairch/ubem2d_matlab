if (~exist('plots','dir'))
    mkdir('plots');
end
% Plot body
data = load('data/cylinder_body.txt');
x = data(:,1);
y = data(:,2);
npan = length(x)-1;
figure;
hold on;
fill(x,y,[.9 .9 .9]);
plot(x,y,'-ok','LineWidth',1.5);
axis('equal');
grid on;
box on;
xlabel('$x$','Interpreter','latex');
ylabel('$y$','Interpreter','latex');
title(sprintf('Cylinder, %d panels',npan),'Interpreter','latex');
set(gca,'FontSize',14);
saveTightFigure('plots/cylinder_body.pdf');

% Plot pressure distribution
data = sortrows(load('data/cylinder_cp.txt'));
th = data(:,1);
Cp = data(:,2);
figure;
plot(th,Cp,'ok',th,1-4*sin(th).^2,'-k','LineWidth',1.5);
axis('equal');
grid on;
box on;
xlabel('$\theta$','Interpreter','latex');
ylabel('$C_p$','Interpreter','latex');
title(['Pressure distribution for uniform flow past a cylinder',...
    sprintf(', %d panels',npan)],'Interpreter','latex');
legend('Numerical','Analytical');
set(gca,'FontSize',14);
saveTightFigure('plots/cylinder_cp.pdf');