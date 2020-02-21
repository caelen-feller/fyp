clf
ff = u{3};
sx = u{1}; sy = u{2};
[xx,yy] = meshgrid(sx,sy); zz = xx + 1i*yy;
true_sol = real(exp(zz));


% Error Heatmaps, both with absolute and relative error
bound = max(max(abs(true_sol)));
imagesc(flip((ff - true_sol)./abs(true_sol))), colorbar
imagesc(flip((ff - true_sol))), colorbar
colormap hot

% Angle plot
per_w = w([1:end 1]);
w_subdiv = [];
for i = 1:length(w)
w_subdiv = [w_subdiv linspace(per_w(i), per_w(i+1),10)];
end
w_subdiv = (w_subdiv - w_c).*0.995 + w_c;

clf
imagesc([0.01,1.99],[1.99,0.01], flipud((ff-true_sol)*10^11./(abs(true_sol) +1)))
set(gca,'YDir','normal')
hold on
x = real(w_subdiv);
y = imag(w_subdiv);
z = zeros(size(x));
col = angle((x+1i*y)-w_c);  % This is the color, vary with x in this case.
surface([x;x],[y;y],[z;z],[col;col],...
        'facecol','no',...
        'edgecol','interp',...
        'linew',5), colorbar;
hold on;
plot(real(w_c), imag(w_c), '.k', "markersize",20)
caxis([-pi pi])


clf
imagesc([0,2],[2,0],flipud(angle(zz-w_c))), colorbar('Ticks',[-0.995*pi, -0.5*pi,0, 0.5*pi,0.995*pi],...
         'TickLabels',{'$-\pi$','$-1/2\pi$','$0$','$1/2\pi$','$\pi$'} ,'TickLabelInterpreter','latex')
set(gca,'YDir','normal')
hold on
plot(w([1:end 1]), "linewidth", 2, 'color', [1 0 (angle(zz-w_c)+pi)/max((angle(zz-w_c)+pi))])


% Overlay code
ax2 = axes;
contour(ax2,ff)
linkaxes([ax1,ax2])
ax2.Visibless='off';
ax2.XTick=[]
ax2.XTick=[];
ax2.YTick=[];
colormap(ax2,'cool')
colormap(ax1,'hot')
set([ax1,ax2],'Position',[.17 .11 .685 .815]);
cb1 = colorbar(ax1,'Position',[.05 .11 .0675 .815]);
cb2 = colorbar(ax2,'Position',[.88 .11 .0675 .815]);


% 