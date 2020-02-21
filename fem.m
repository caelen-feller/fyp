model = createpde();
L = [2 6 0 2 2 1 1 0 0 0 1 1 2 2];
dl = decsg(L.');
geometryFromEdges(model, dl);
pdegplot(model,'EdgeLabels','on'); 
axis equal

clf
bc = @(location,state) real(exp(location.x + 1i*location.y));
applyBoundaryCondition(model,'dirichlet','Edge',1:model.Geometry.NumEdges,'r', bc);
specifyCoefficients(model,'m',0,'d',0,'c',-1,'a',0,'f',0);
hmax = 0.1;
generateMesh(model,'Hmax',hmax);

figsize = [0,0,2.5,2.5];
fontsize = 9;
h = figure;
pdemesh(model)
xlim([0 2]); ylim([0 2]);
set(h, 'Units','Inches');
pos = get(h,'Position');
set(h,...
    'Position', figsize,...
    'Units', 'Inches', ...
    'PaperPositionMode','Auto',...
    'PaperUnits','Inches')
set(gca, 'Units', 'normalized',...
    'FontUnits', 'points',...
    'FontWeight', 'normal',...
    'FontSize', fontsize, ...
    'FontName', 'Times',...
    'TickLabelInterpreter', 'latex');
print(h,'report\figures\fem-mesh','-depsc')


results = solvepde(model);
u = results.NodalSolution;
pdeplot(model,'XYData',u)
title('Numerical Solution');
xlabel('x')
ylabel('y')

figsize = [0,0,3,2.5];
fontsize = 9;
h = figure;
p = model.Mesh.Nodes;
exact = real(exp(p(1,:) + 1i.*p(2,:)));
pdeplot(model, 'XYData', u - exact.');

xlim([0 2]); ylim([0 2]);
set(h, 'Units','Inches');
pos = get(h,'Position');
set(h,...
    'Position', figsize,...
    'Units', 'Inches', ...
    'PaperPositionMode','Auto',...
    'PaperUnits','Inches')
set(gca, 'Units', 'normalized',...
    'FontUnits', 'points',...
    'FontWeight', 'normal',...
    'FontSize', fontsize, ...
    'FontName', 'Times',...
    'TickLabelInterpreter', 'latex');
print(h,'report\figures\fem-error','-depsc')

hmax = 0.1;
errs = [];
err = 1;
while err > 1e-8 % run until error <= 5e-7
    generateMesh(model,'Hmax',hmax); % refine mesh
    results = solvepde(model);
    u = results.NodalSolution;
    p = model.Mesh.Nodes;
    exact = real(exp(p(1,:) + 1i.*p(2,:)));
    err = norm(u - exact',inf); % compare with exact solution
    disp("err = " + err + " hmax = " +hmax);
    
    errs = [errs err]; % keep history of err
    hmax = hmax/2;
end

semilogy(  sqrt(1:numel(errs)), errs,'MarkerSize',12);
ax = gca;
ax.XTick = 1:numel(errs);
title('Error History');
xlabel('\sqrt{order}');
ylabel('Norm of Error');

p = model.Mesh.Nodes;
exact = (1 - p(1,:).^2 - p(2,:).^2)/4;
pdeplot(model,'XYData',u - exact')
title('Error');
xlabel('x')
ylabel('y')