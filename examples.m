% Example of use of framework on various domains
% Define example polygonal domains
w_square = [1-1i; 1+1i; 1i-1; -1-1i];
w_L = [2; 2+1i; 1+1i; 1+2i; 2i; 0];
w_isodrum = ([1+2i; 1+3i; 2i; 1i+1; 2+1i; 2; 3+1i; 3+2i]-(1.5+1.5i))/1.8;
w_long_no_fix = [-4 4 4+1i .5i -4+1i].';
w_long = [linspace(-4,4,10) (4+1i + (.5i - 4-1i).*linspace(0,1,10)) (0.5i + (-4+1i - .5i).*linspace(0.1,1,9))].';
w_spike = [0 1 1+.3i .5+.5i 1+.7i 1+1i 1i].';
w_rand = (exp(2i*pi*(1:10)/10).*(.1+rand(1,10))).';

% This domain
% plot(exp(1i*(linspace(0,1.5*pi,100) + pi/4))), hold on
% plot(exp(1i*(5*pi/4 - linspace(0,pi/2,100))) + sqrt(2)), hold off
% Give corners and arclength param of edges for curved domain
w_circular = [exp(1.75i*pi); exp(.25i*pi)];
circular = [{@(t) exp(1i*(t+pi/4))} 1.5*pi 1; ...
    {@(t) exp(1i*(5*pi/4 - t)) + sqrt(2)} pi/2 1];

w_tri = [exp(.25i*pi); -1; exp(1.75i*pi); -0.5];
tri = [{@(t) exp(1i*(t+pi/4))} 1.5*pi/2 1; ...
    {@(t) exp(1i*(t+pi/4+1.5*pi/2))} 1.5*pi/2 -1; ...
    {@(t) exp(1.75i*pi)- t*(0.5 + exp(1.75i*pi))} 1 1; ...
    {@(t) -0.5 + t*(exp(.25i*pi) + 0.5)} 1 1];

curved = false; curves = [];
w = w_L;

% Need boundary condition, specified as real funcs of boundary values
% Other possible ways of specifiying this data include:
%  - Real funcs on arc length
%  - Array of constant values for the sides

h_linear = repmat({@(z) real(z)}, 1, length(w));
h_quad = repmat({@(z) real(z)^2}, 1, length(w));
h_abs = repmat({@(z) abs(z)}, 1, length(w));
h_sin = repmat({@(z) sin(4*real(z))}, 1, length(w));
h_discont = [repmat({@(z) 0}, 1, length(w)/2) ...
   repmat({@(z) 1}, 1, length(w)/2);];
h_robin = [repmat([{@(z) abs(z)} 1; {@(z) abs(z)} 0], length(w)/2, 1)];
h = h_robin; weighted = true; neumann = true;

curved_hull = [];
if curved
    [m, ~] = size(w);
    for i = 1:m
        prev = mod(m + i - 2, m)+1; next = mod(i, m)+1;
        curved_hull = [ curved_hull ...
            curves{i,1}(linspace(0,curves{i,2},100))];
    end
    curved_hull = curved_hull.';
end

% Center of domain (used as sample point for the polynomial terms)
w_c = 0;

tol = 1e-4 % Desired accuracy of solution

u = laplace_solver(w, w_c, h, tol...
    ,'tests', false ...
    ,'neumann', neumann ...
    ,'plot3', false ...
    ,'discont', weighted ...
    ,'curved', curved ...
    ,'curves', curves ...
    ,'curved_hull', curved_hull);
