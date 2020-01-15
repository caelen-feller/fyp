% Example of use of framework on the L-Shaped Domain

% Define example polygonal domains
w_square = [1-1i; 1+1i; 1i-1; -1-1i];
w_L = [2; 2+1i; 1+1i; 1+2i; 2i; 0];
w_isodrum = ([1+2i; 1+3i; 2i; 1i+1; 2+1i; 2; 3+1i; 3+2i]-(1.5+1.5i))/1.8; 

% This domain
% plot(exp(1i*(linspace(0,1.5*pi,100) + pi/4))), hold on
% plot(exp(1i*(5*pi/4 - linspace(0,pi/2,100))) + sqrt(2)), hold off
% Give corners and arclength param of edges for curved domain
curved = false; curves = [];
w_curved = [exp(1.75i*pi); exp(.25i*pi)];
curves = [{@(t) exp(1i*(t+pi/4))} 1.5*pi 1; ...
    {@(t) exp(1i*(5*pi/4 - t)) + sqrt(2)} pi/2 1];
w = w_L;

% Need boundary condition, specified as real funcs of boundary values
% Other possible ways of specifiying this data include:
%  - Real funcs on arc length
%  - Array of constant values for the sides

h_linear = repmat({@(z) real(z)}, 1, length(w));
h_quad = repmat({@(z) real(z)^2}, 1, length(w));
h_sin = repmat({@(z) sin(5*real(z))}, 1, length(w));
h_discont = [repmat({@(z) 0}, 1, length(w)/2) ...
    repmat({@(z) 1}, 1, length(w)/2);];
h = h_sin; weighted = false;

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
w_c = 0.5 + 0.5i;
if curved, w_c = mean(curved_hull); end
    
tol = 1e-6 % Desired accuracy of solution 

u = laplace_solver(w, w_c, h, tol...
    ,'tests', false ...
    ,'plot3', false ...
    ,'discont', weighted ...
    ,'curved', curved ...
    ,'curves', curves ...
    ,'curved_hull', curved_hull);

