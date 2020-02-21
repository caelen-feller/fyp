% Example of use of framework on various domains
% Define example polygonal domains
% Note: not all bcs will work on domains with non-even number of corners
w_square = [1-1i; 1+1i; 1i-1; -1-1i];
w_L = [2; 2+1i; 1+1i; 1+2i; 2i; 0];
w_isodrum = ([1+2i; 1+3i; 2i; 1i+1; 2+1i; 2; 3+1i; 3+2i]-(1.5+1.5i))/1.8;
w_long_no_fix = [-4 4 4+1i .5i -4+1i].';
w_long = [linspace(-4,4,10) (4+1i + (.5i - 4-1i).*linspace(0,1,9)) (0.5i + (-4+1i - .5i).*linspace(0.1,1,9))].';

mag= 0.75;
w_spike = [0 1 1+ 0.3i (1-mag)*1+.5i 1+ .7i 1+1i 1i 0.5i].';

w_slit = [0 1 1+ mag*0.5i .5+.5i 1+ (2-mag)*.5i 1+1i 1i 0.5i].';

w_rev_spike = [0 1 1+.3i 1.5+.5i 1+.7i 1+1i 1i].';
w_rand = (exp(2i*pi*(1:10)/10).*(.1+rand(1,10))).';

w_star_a = (exp(2i*pi*(1:10).'/10)).*(.5*exp(2i*pi/20));
w_star_b = (exp(2i*pi*(1:10).'/10));
w_star = reshape([w_star_a(:) w_star_b(:)]',2*size(w_star_a,1), [])

w_star_2_a = (exp(2i*pi*(1:10).'/10)).*(exp(2i*pi*15/100));
w_star_2_b = (exp(2i*pi*(1:10).'/10));
w_star_2 = reshape([w_star_2_a(:) w_star_2_b(:)]',2*size(w_star_a,1), [])

% Give corners and arclength param of edges for curved domains
w_circular = [exp(1.75i*pi); -1; exp(.25i*pi); sqrt(2)-1];
circular = [{@(t) exp(1i*(t+pi/4))} 1.5*pi/2 1; ...
    {@(t) exp(1i*(t+pi/4+1.5*pi/2))} 1.5*pi/2 -1; ...
    {@(t) exp(1i*(5*pi/4 - t)) + sqrt(2)} pi/4 1; ... 
    {@(t) exp(1i*(pi - t)) + sqrt(2)} pi/4 1;];

w_tri = [exp(.25i*pi); -1; exp(1.75i*pi); -0.5];
tri = [{@(t) exp(1i*(t+pi/4))} 1.5*pi/2 1; ...
    {@(t) exp(1i*(t+pi/4+1.5*pi/2))} 1.5*pi/2 -1; ...
    {@(t) exp(1.75i*pi)- t*(0.5 + exp(1.75i*pi))} 1 1; ...
    {@(t) -0.5 + t*(exp(.25i*pi) + 0.5)} 1 1];

curved = true; curves = circular;
% curved = true; curves = tri; % example for it on with w_tri
w = w_circular;

% Need boundary conditions, specified as real funcs of boundary values
% Could also specify as set of functions on arc length

h_linear = repmat({@(z) real(z)}, length(w), 1);
h_quad = repmat({@(z) real(z)^2}, length(w), 1);
h_abs = repmat({@(z) abs(z)}, length(w), 1);
h_sin = repmat({@(z) real(exp(z))}, length(w), 1);
h_blowup = [ {@(z) 1/(imag(z)-1.1)}; {@(z) 1/(real(z)-1.1)} ;{@(z)1/(real(z)-1.1)} ;{@(z) 1/(imag(z)-1.1)}];

% These all require weighted=true to converge properly
h_discont = [repmat({@(z) 0}, length(w)/2, 1);...
    repmat({@(z) 1}, length(w)/2, 1);];
% % These require neumann support to be on
h_r_simple = [{@(z) 1} 1; repmat([{@(z) 0} 0], length(w) - 1, 1)];
h_r_demo = [{@(z) 1} 1; {@(z) 0} 0; {@(z) 0} 0; {@(z) 1} 1; ...
    {@(z) 0} 0; {@(z) 0} 0;];


h_r_test = [{@(z) 1} 1; {@(z) real(z)} 0; {@(z) -1} 0; {@(z) real(z)} 0];
weighted = true; neumann = false;
h = h_sin; 

% Creates a uniformly sampled hull of the boundary for curved domains
% Used in inclusion testing and plotting
n_slices = 100;
curved_hull = [];
if curved
    [m, ~] = size(w);
    for i = 1:m
        prev = mod(m + i - 2, m)+1; next = mod(i, m)+1;
        curved_hull = [ curved_hull ...
            curves{i,1}(linspace(0,curves{i,2},n_slices))];
    end
    curved_hull = curved_hull.';
end

% Center of domain (used as sample point for the polynomial terms)
w_c =0;

tol = 1e-16 % Desired accuracy of solution

u = laplace_solver(w, w_c, h, tol...
    ,'tests', false ...
    ,'neumann', neumann ...
    ,'plot3', false ...
    ,'heatmap', true ...
    ,'solfunc', false ...
    ,'discont', weighted ...
    ,'curved', curved ...
    ,'curves', curves ...
    ,'curved_hull', curved_hull);
