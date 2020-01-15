% Laplace equation solver on polygonal domains
% Arguments:
%   w = List of boundary vertices, counterclockwise
%   center = Point centrally situated within domain used for interior
%   polynomial-based approximation
%   h = Boundary value problem (Dirichlet), specified per-edge as functions
%   of the points on the boundary
%   tol = Desired accuracy of solution
% Return Value: 
%   u = Solution to boundary value problem, valid on the domain only
function u = laplace_solver(w, center, h, tol, varargin)
    p = inputParser;
    addOptional(p, 'tests', false);
    addOptional(p, 'neumann', false);
    addOptional(p, 'discont', false);
    addOptional(p, 'curved', false)
    addOptional(p, 'curves', []);
    addOptional(p, 'curved_hull', []);
    addOptional(p, 'plot3', false);
    
    parse(p, varargin{:});
    curved = p.Results.curved;
    if curved, curves = p.Results.curves; end

    % Get the maximal diameter of the domain
    w_r = sort(real(w)); w_r = w_r([1 end]);
    w_i = sort(imag(w)); w_i = w_i([1 end]);
    scale = max([diff(w_r),diff(w_i)]);
    
    if curved
        b_r = sort(real(p.Results.curved_hull)); b_r = b_r([1 end]);
        b_i = sort(imag(p.Results.curved_hull)); b_i = b_i([1 end]);
        scale = max([diff(b_r),diff(b_i)]);
    end
        
    % Confirm validity of poles and boundary sampling points if tests on
    if p.Results.tests
        poles = compute_poles_test(w, scale, 3, center);
        bdd = bdd_sample_test(w, poles, h, 3);
    end
    
    % Increase n until convergence to specified tolerance.
    prev_err = inf; bdd = []; 
    for n = 3:50
        % Exp. clustered poles at corners  
        poles = compute_poles(w, scale, n, p);
        flat_poles = cell2mat(poles).';
        [~, N1] = size(flat_poles);
        % Distance of poles from corners and boundary sample points
        [bdd, dist] = bdd_sample(w, poles, h, n, p); 
        [M,~] = size(bdd); 
        b = bdd(:,2); % Values at boundary points
        newman = dist./(bdd(:,1) - flat_poles);

        % Polynomial approx in bulk
        N2 = ceil(n/2);
        runge = ((bdd(:,1) - center)./scale).^(1:N2); 
        
        % Construct least squares matrix
        N = 2*N1 + 2*N2;
        disp("n = "+n); disp("N = "+N); disp("N1 = "+N1); disp("N2 = "+N2); disp("M = "+M);
        A = [real(runge) imag(runge) real(newman) imag(newman) ones(M,1)];
        
        % If necessary, weight by distance to vertex
        weights = ones(1,M);
        if p.Results.discont, weights = min(abs(bdd(:,1).' - w)); end
        weight_m = spdiags(sqrt(weights.'), 0, M, M);
        
        % Solve using backslash, matricies are always ill-conditioned
        warn = warning('off','MATLAB:rankDeficientMatrix');
        c = (weight_m * A) \ (weight_m * b);
        warning(warn.state,'MATLAB:rankDeficientMatrix');
        
        % Check error on boundary
        err = A*c - b;
        err = norm(weights.*err,inf);
        disp("err="+err);
                
        % Exit if tolerance reached or if error begins to increase 
        % as convergence is not possible due to numerical instability.
        if err <= tol || err > prev_err*10 
            break;
        end
        prev_err = err;
    end
    
    % Construct proposed solution
    u = @(z) [real(((z-center)./scale).^(1:N2)) ...
            imag(((z-center)./scale).^(1:N2))   ... 
            real(dist./(z-flat_poles))        ...
            imag(dist./(z-flat_poles)) ones(size(z))]*c;
    
    % Plot solution, poles and boundary sampling points   
    % Plotting Config
    LW = 'linewidth'; MS = 'markersize'; FS = 'fontsize';
    fs = 8; PO = 'position'; FW = 'fontweight'; NO = 'normal';
    
    if curved, w_r = b_r; w_i = b_i; w = p.Results.curved_hull; end
    ax = [w_r(1:2); w_i(1:2)] + .2*scale*[-1 1 -1 1]';
    axwide = [w_r(1:2); w_i(1:2)] + 1.1*scale*[-1 1 -1 1]';

    % Create sampling grid
    sx = linspace(w_r(1),w_r(2),100); sy = linspace(w_i(1),w_i(2),100);
    [xx,yy] = meshgrid(sx,sy); zz = xx + 1i*yy;
    
    % Evaluate solution on grid, discard points not in domain
    ff = zz; ff(:) = u(zz(:)); ff(~indomain(zz,w,p)) = NaN;
    
    % Get bounds for color scale, plot solution, domain and poles/bdd
    clf, shg
    axes(PO,[.4 .25 .6 .6])
    levels = linspace(min(min(real(ff))),max(max(real(ff))),20);
    if p.Results.plot3, plot3(xx,yy,ff), colorbar; 
    else
        contour(sx,sy,ff,levels,LW,.5), colorbar, axis equal, hold on;
        plot(w([1:end 1]), '-k', LW,1), plot(flat_poles, '.r', MS,5);
        set(gca,FS,fs-1), plot(real(center),imag(center),'.k',MS,6), axis(ax)
    end
    title(['#bdd sample points = ' int2str(M) ...
        ', #poles = ' int2str(N1)],FS,fs,FW,NO), hold off;
    
    % Plot boundary error
    axes(PO,[.05 .4 .25 .2])
    errmin = .01*tol; ws = 'error'; if p.Results.discont, ws = 'weighted error'; end
    axis([-pi pi .0001*errmin 1]), grid on
    semilogy(angle(bdd(:,1)-center),weights.*abs(u(bdd(:,1))-bdd(:,2)),'.k',MS,4), hold off
    set(gca,'ytick',10.^(-16:4:0))
    set(gca,'xtick',pi*(-1:1),'xticklabel',{'-\pi','0','\pi'})
    set(gca,FS,fs-1), xlabel('angle on boundary wrt wc',FS,fs)
    title([ws ' on boundary'],FS,fs,FW,NO)

    % TODO: Plot speed of convergence 
    % TODO: A posteriori error on finer boundary sampling
    
end

% Utility to check if point is in domain
function in = indomain(z, w, p)
    % TODO: make this accurate for curved boundaries
    if p.Results.curved, w = p.Results.curved_hull; end
    in = inpolygon(real(z),imag(z),real(w),imag(w)); 
end 

% Generate exponentially clustered poles along vertex normals of domain
% TODO: Vectorize computations further
function poles = compute_poles(w, scale, n, p)
    [m, ~] = size(w);
    curved = p.Results.curved;
    if curved, curves = p.Results.curves; end

    % Edge normals
    n_e = zeros(m, 1);
    
    for i = 1:m
        n_e(i) = (w(i) - w(mod(i + m - 2, m)+1))*-1i;
        n_e(i) = n_e(i)/norm(n_e(i));
    end
    
    % Vertex normals (external angle bisectors)
    n_v = zeros(m, 1);

    for i = 1:m
        n_v(i) = n_e(i) + n_e(mod(i, m)+1);
        n_v(i) = n_v(i)/norm(n_v(i));
    end
    
    % Pole generation 
    poles = cell(m,1);
    sigma = 4;
    re_entrant = real(n_v).*real(n_e*1i)+imag(n_v).*imag(n_e*1i);
    if curved
        re_entrant = -ones(m,1); 
        n_v = cell2mat(curves(:,3));
    end
    for i = 1:m
        N = n;
        % Check for re-entrancy 
        if re_entrant(i) < 0
            N = 3*n;
        end
        
        dist = exp(-sigma*(sqrt(N) - sqrt(1:N)))';
        pole = w(i) + n_v(i) * scale .* dist;
        % Check for poles in domain
        poles{i} = pole(~indomain(pole,w,p));
    end
end

% Generate boundary sample points and get pole scale (distance to corners)
% TODO: Vectorize computations further
function [bdd,dist] = bdd_sample(w, poles, h, n, p)
    [m, ~] = size(w);
    bdd = zeros(3*n*m,2);
    dist = zeros(3*n*m,1);
    bdd_i=1; dist_i=1;
    curved = p.Results.curved;
    if curved, curves = p.Results.curves; end
    
    for i = 1:m
        [N, ~] = size(poles{i}); 
        prev = mod(m + i - 2, m)+1; next = mod(i, m)+1;
        neg = w(prev) - w(i); pos = w(next) - w(i);
        if curved, neg = curves{prev,2}; pos = curves{i,2}; end
        % loop through poles, checking distance
        for j = 1:N
            dist(dist_i) = norm(w(i) - poles{i}(j)); 
            for k = 1:3
                t = k/3 * dist(dist_i);
                if(t < norm(neg))
                    bdd(bdd_i, 1) = t * neg/norm(neg) + w(i);
                    if curved, bdd(bdd_i, 1) = ...
                        curves{prev, 1}(curves{prev, 2} - t); end
                    bdd(bdd_i, 2) = h{prev}(bdd(bdd_i,1));
                    bdd_i = bdd_i + 1;
                end
                if(t < norm(pos))
                    bdd(bdd_i, 1) = t * pos/norm(pos) + w(i);
                    if curved, bdd(bdd_i, 1) = curves{i, 1}(t); end
                    bdd(bdd_i, 2) = h{i}(bdd(bdd_i,1));
                    bdd_i = bdd_i + 1;
                end
            end
            dist_i = dist_i+1;
        end
    end
    bdd = bdd(1:(bdd_i-1),:);
    dist = dist(1:(dist_i-1),:).';
end

% Tests for the above functions
function out = compute_poles_test(w, scale, n, center)
    LW = 'linewidth'; MS = 'markersize'; 
    plot(w([1:end 1]), '-k', LW,1); 
    out = compute_poles(w, scale, n);
    flat_poles = cell2mat(out);
    hold on; plot(flat_poles, '.g', MS,5); hold off
end

function out = bdd_sample_test(w, poles, h, n)
    out = bdd_sample(w, poles, h, n);
    LW = 'linewidth'; MS = 'markersize';  
    hold on; plot(out(:,1), '.g', MS,5); hold off
    disp(' ')
    y = input('next? ');
end
