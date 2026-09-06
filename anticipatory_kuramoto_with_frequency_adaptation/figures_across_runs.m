printfigureflag = 0;

%% A)
G = groupsummary(REZ, {'k2','taus2'}, 'mean', 'deltaOmega');

k2_dW    = unique(G.k2);
taus2_dW = unique(G.taus2);
[K2_dW, TAUS2_dW] = meshgrid(k2_dW, taus2_dW);

Z_dW = nan(size(K2_dW));
for i = 1:height(G)
    r = find(taus2_dW == G.taus2(i));
    c = find(k2_dW    == G.k2(i));
    Z_dW(r,c) = G.mean_deltaOmega(i);
end

figure(1)
surf(TAUS2_dW, K2_dW, Z_dW, 'EdgeColor', 'none');
colormap('cool');
cb = colorbar;
set(gca, 'XDir', 'reverse', 'YDir', 'reverse');
hold on
Z_zero_dW = zeros(size(K2_dW));
m = mesh(TAUS2_dW, K2_dW, Z_zero_dW, 'EdgeColor', 'k', 'FaceColor', 'none', 'FaceAlpha', 0);
contour3(TAUS2_dW, K2_dW, Z_dW, [0 0], 'k', 'LineWidth', 2);
hold off
xlabel('\tau_{follower,2} [degrees]')
ylabel('k_2')
ylabel(cb, '\Delta \omega [%]')
xlim([0 2*pi])
ylim([0 4])
set(gca,'XTick',round(0:pi/2:2*pi,2))
set(gca,'XTickLabel',(0:pi/2:2*pi)./pi*180)
zlabel('\Delta \omega, %')
set(gcf,'color','w')
set(gcf,'PaperPosition',[0 0 5 4])
set(gcf,'InvertHardcopy','off')
set(gca,'fontsize',14)
if printfigureflag == 1
    set(gca,'view',[156 20])
    f = fullfile(pwd,['omega_drift_mean_3d_' char(datetime('now','TimeZone','local','Format','y-MM-d-hhmmss'))]);
    print('-djpeg','-r600',[f '.jpg'])
    print('-dsvg',[f '.svg'])

    set(gca,'view',[180 90])
    f = fullfile(pwd,['omega_drift_2d_' char(datetime('now','TimeZone','local','Format','y-MM-d-hhmmss'))]);
    print('-djpeg','-r600',[f '.jpg'])
    print('-dsvg',[f '.svg'])
end



%% B)
G = groupsummary(REZ, {'k2','taus2'}, 'mean', 'CV');

k2_dW    = unique(G.k2);
taus2_dW = unique(G.taus2);
[K2_dW, TAUS2_dW] = meshgrid(k2_dW, taus2_dW);

Z_dW = nan(size(K2_dW));
for i = 1:height(G)
    r = find(taus2_dW == G.taus2(i));
    c = find(k2_dW    == G.k2(i));
    Z_dW(r,c) = G.mean_CV(i);
end
baselinecv = mean(Z_dW(K2_dW==0));

figure(2)
surf(TAUS2_dW, K2_dW, Z_dW, 'EdgeColor', 'none');
colormap('cool');
cb = colorbar;
set(gca, 'XDir', 'reverse', 'YDir', 'reverse');
hold on
Z_zero_dW = zeros(size(K2_dW)) + baselinecv;
m = mesh(TAUS2_dW, K2_dW, Z_zero_dW, 'EdgeColor', 'k', 'FaceColor', 'none', 'FaceAlpha', 0);
contour3(TAUS2_dW, K2_dW, Z_dW, [baselinecv baselinecv], 'k', 'LineWidth', 2);
hold off
xlabel('\tau_{follower,2} [degrees]')
ylabel('k_2')
ylabel(cb, 'CV_{IOI}, %')
xlim([0 2*pi])
ylim([0 4])
set(gca,'XTick',round(0:pi/2:2*pi,2))
set(gca,'XTickLabel',(0:pi/2:2*pi)./pi*180)
zlabel('CV_{IOI}, %')
set(gcf,'color','w')
set(gcf,'PaperPosition',[0 0 5 4])
set(gcf,'InvertHardcopy','off')
set(gca,'fontsize',14)
if printfigureflag == 1
    set(gca,'view',[157 11])
    f = fullfile(pwd,['cv_mean_3d_' char(datetime('now','TimeZone','local','Format','y-MM-d-hhmmss'))]);
    print('-djpeg','-r600',[f '.jpg'])
    print('-dsvg',[f '.svg'])

    set(gca,'view',[180 90])
    f = fullfile(pwd,['cv_mean_2d_' char(datetime('now','TimeZone','local','Format','y-MM-d-hhmmss'))]);
    print('-djpeg','-r600',[f '.jpg'])
    print('-dsvg',[f '.svg'])
end


