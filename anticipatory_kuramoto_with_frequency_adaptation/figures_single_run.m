printfigureflag = 0;

figure(1)
surf(X,Y,reshape(cycle_variability*100,size(X,1),[]),"EdgeColor","black")
xlabel('\tau_2 [rad]')
ylabel('k_2')
zlabel('CV_{IOI}, %')
xlim([0 2*pi])
ylim([0 4])
set(gcf,'color','w')
set(gcf, 'PaperPosition', [0 0 4 3])
set(gcf, 'InvertHardcopy', 'off')
set(gca,'fontsize',14)
if printfigureflag == 1
    f = fullfile(pwd,['rel_phase_final_' char(datetime('now','TimeZone','local','Format','y-MM-d-hhmmss')) '.jpeg']);
    print('-djpeg','-r600',f)
end


figure(2)
surf(X,Y,reshape(cycle_variability*100,size(X,1),[]),"EdgeColor","none")
colormap cool
cb = colorbar;
ylabel(cb, 'CV_{IOI}, %')
xlabel('\tau_{follower,2} [degrees]')
ylabel('k_2')
xlim([0 2*pi])
ylim([0 4])
set(gca, 'XDir', 'reverse')
set(gca, 'YDir', 'reverse')
set(gca,'XTick',round(0:pi/2:2*pi,2))
set(gca,'XTickLabel',(0:pi/2:2*pi)./pi*180)
set(gca,'YTick',round(0:1:4))
set(gca,'YTickLabel',round(0:1:4))
set(gca,'view',[180 90]) % set(gca,'view',[-266 5])
set(gcf,'color','w')
set(gcf, 'PaperPosition', [0 0 5 4])
set(gcf, 'InvertHardcopy', 'off')
set(gca,'fontsize',14)
if printfigureflag == 1
    f = fullfile(pwd,['omega_drift_' char(datetime('now','TimeZone','local','Format','y-MM-d-hhmmss')) '.jpeg']);
    print('-djpeg','-r600',f)
end


figure(3)
surf(X,Y,reshape(omega_delta,size(X,1),[]),"EdgeColor","none")
colormap cool
cb = colorbar;
ylabel(cb, '\Delta \omega [%]')
xlabel('\tau_{follower,2} [degrees]')
ylabel('k_2')
xlim([0 2*pi])
ylim([0 4])
set(gca, 'XDir', 'reverse')
set(gca, 'YDir', 'reverse')
set(gca,'XTick',round(0:pi/2:2*pi,2))
set(gca,'XTickLabel',(0:pi/2:2*pi)./pi*180)
set(gca,'YTick',round(0:1:4))
set(gca,'YTickLabel',round(0:1:4))
set(gca,'view',[180 90]) % set(gca,'view',[-266 5])
zlabel('\Delta \omega, %')
set(gcf,'color','w')
set(gcf, 'PaperPosition', [0 0 5 4])
set(gcf, 'InvertHardcopy', 'off')
set(gca,'fontsize',14)
if printfigureflag == 1
    f = fullfile(pwd,['omega_drift_' char(datetime('now','TimeZone','local','Format','y-MM-d-hhmmss')) '.jpeg']);
    print('-djpeg','-r600',f)
end


figure(4)
index = k_follower_2_vec == k_follower_1;
% tau_diff = (tau_follower_1 + tau_follower_2_vec)./2 - tau_leader;
tau_diff = tau_follower_2_vec;
plot(tau_diff(index),omega_delta(index),'-ok','linewidth',2);
hold on
plot(tau_diff(index),tau_diff(index)*0,'--k','linewidth',1);
hold off
ylabel('\Delta \omega [%]')
% xlabel('\tau_{follower} - \tau_{leader} [degrees]')
xlabel('\tau_{follower} [degrees]')
xlim([min(min(0,tau_diff)) max(tau_diff)])
set(gca,'XTick',round(0:pi/2:2*pi,2))
set(gca,'XTickLabel',(0:pi/2:2*pi)./2/pi*360)
set(gcf,'color','w')
set(gcf, 'PaperPosition', [0 0 5 4])
set(gcf, 'InvertHardcopy', 'off')
set(gca,'fontsize',14)
if printfigureflag == 1
    f = fullfile(pwd,['omega_drift_by_tau_' char(datetime('now','TimeZone','local','Format','y-MM-d-hhmmss')) '.jpeg']);
    print('-djpeg','-r600',f)
end


figure(5)
[~,index] = min(abs(tau_follower_2_vec - pi/2));
baseline_tau2 = tau_follower_2_vec(index);
index = tau_follower_2_vec == baseline_tau2;
plot(k_follower_2_vec(index,1),omega_delta(index,1),'-ok','linewidth',2);
hold on
plot(k_follower_2_vec(index,1),k_follower_2_vec(index,1)*0,'--k','linewidth',1);
hold off
ylabel('\Delta \omega [%]')
xlabel('k_2')
xlim([0 4])
set(gca,'XTick',0:4)
set(gca,'XTickLabel',0:4)
set(gcf,'color','w')
set(gcf, 'PaperPosition', [0 0 5 4])
set(gcf, 'InvertHardcopy', 'off')
set(gca,'fontsize',14)
if printfigureflag == 1
    f = fullfile(pwd,['omega_drift_by_phi0_' char(datetime('now','TimeZone','local','Format','y-MM-d-hhmmss')) '.jpeg']);
    print('-djpeg','-r600',f)
end


figure(6)
[~,index] = min(abs(tau_follower_2_vec - tau_follower_1));
baseline_tau1 = tau_follower_2_vec(index);
index = tau_follower_2_vec == baseline_tau1;
plot(k_follower_2_vec(index,1),cycle_variability(index,1)*1e2,'-ok','linewidth',2);
hold on
plot(k_follower_2_vec(index,1),...
    k_follower_2_vec(index,1)*0 + cycle_variability(k_follower_2_vec==0 & tau_follower_2_vec == baseline_tau1,1)*1e2,...
    '--k','linewidth',1);
hold off
ylabel('CV_{IOI,follower}, %')
xlabel('k_2')
xlim([0 4])
set(gca,'XTick',0:4)
set(gca,'XTickLabel',0:4)
set(gcf,'color','w')
set(gcf, 'PaperPosition', [0 0 5 4])
set(gcf, 'InvertHardcopy', 'off')
set(gca,'fontsize',14)
if printfigureflag == 1
    f = fullfile(pwd,['omega_drift_by_phi0_' char(datetime('now','TimeZone','local','Format','y-MM-d-hhmmss')) '.jpeg']);
    print('-djpeg','-r600',f)
end


figure(7)
[~,index] = min(abs(tau_follower_2_vec - tau_follower_1));
baseline_tau1 = tau_follower_2_vec(index);
index = k_follower_2_vec == k_follower_1;
tau_diff = tau_follower_2_vec;
plot(tau_diff(index),cycle_variability(index)*1e2,'-ok','linewidth',2);
hold on
plot(tau_diff(index,1),...
    tau_diff(index,1)*0 + cycle_variability(k_follower_2_vec==0 & tau_follower_2_vec == baseline_tau1,1)*1e2,...
    '--k','linewidth',1);
hold off
ylabel('CV_{IOI} [%]')
xlabel('\tau_{follower} [degrees]')
xlim([min(min(0,tau_diff)) max(tau_diff)])
set(gca,'XTick',round(0:pi/2:2*pi,2))
set(gca,'XTickLabel',(0:pi/2:2*pi)./2/pi*360)
set(gcf,'color','w')
set(gcf, 'PaperPosition', [0 0 5 4])
set(gcf, 'InvertHardcopy', 'off')
set(gca,'fontsize',14)
if printfigureflag == 1
    f = fullfile(pwd,['omega_drift_by_tau_' char(datetime('now','TimeZone','local','Format','y-MM-d-hhmmss')) '.jpeg']);
    print('-djpeg','-r600',f)
end


