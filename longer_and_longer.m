errs = []
for i = 1:10
diameter = i;
w_long_convex = [-diameter-1i diameter-1i diameter+1i -diameter+1i].';
w_long_no_fix = [-diameter diameter diameter+1i .5i -diameter+1i].';
density = i^2;
w_long = [linspace(-diameter,diameter,density) (diameter+1i + (.5i - diameter-1i).*linspace(0,1,density)) (0.5i + (-diameter+1i - .5i).*linspace(0.1,1,density))].';

curved = false; 
w = w_long;

h = repmat({@(z) real(z)^2}, length(w), 1);
weighted = false; neumann = false;

w_c = 0;

tol = 1e-15 % Desired accuracy of solution

u = laplace_solver(w, w_c, h, tol...
    ,'tests', false ...
    , 'plots', false ...
    ,'plot3', false...
    ,'neumann', neumann ...
    ,'solfunc', false ...
    ,'discont', weighted ...
    ,'fixed', true);
errs = [errs u{4}];
end

plot(errs)
