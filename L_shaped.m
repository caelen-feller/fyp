% Example of use of framework on the L-Shaped Domain

% Define example polygonal domains
w_square = [1-1i; 1+1i; 1i-1; -1-1i];
w_L = [2; 2+1i; 1+1i; 1+2i; 2i; 0];
w_isodrum = ([1+2i; 1+3i; 2i; 1i+1; 2+1i; 2; 3+1i; 3+2i]-(1.5+1.5i))/1.8; 

w = w_isodrum;

% Center of domain (used as sample point for the polynomial terms)
w_c = mean(w);

% Need boundary condition, specified as real funcs of boundary values
% Other possible ways of specifiying this data include:
%  - Real funcs on arc length
%  - Array of constant values for the sides

h_linear = repmat({@(z) real(z)}, 1, length(w));
h_quad = repmat({@(z) real(z)^2}, 1, length(w));
h_discont = [repmat({@(z) 0}, 1, length(w)/2) ...
    repmat({@(z) 1}, 1, length(w)/2);];
h = h_discont;

tol = 1e-6 % Desired accuracy of solution 

u = laplace_solver(w, w_c, h, tol, 'discont', true);

