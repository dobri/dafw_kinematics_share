printfigureflag = 0;

% Phase and freq adaptation with follower delays > than leader delay
anticipatory_kuramoto([4 4],[10 10],[pi/5 pi/2],pi/6,0,[],1,printfigureflag);

% Phase adaptation with follower delays > than leader delay
anticipatory_kuramoto([4 0],[10 0],[pi/5 pi/4],pi/6,0,[],1,printfigureflag);

% Parameter space of initial relative phases and tau_2 for fixed K and epsilon
Kvec = [2 2];
epsilon = [5 5];

tau_vec = 0:(pi/10):2*pi;
init_rel_phase = 0:pi/10:2*pi;
[X, Y] = ndgrid(tau_vec,init_rel_phase);
param_space = [X(:) Y(:)];
Z = zeros(size(param_space,1),1);
rel_phase_r = zeros(size(param_space,1),1);
rel_phase = zeros(size(param_space,1),1);
rel_phase_error = zeros(size(param_space,1),1);
tau = zeros(size(param_space,1),2);
tau_1 = 2*pi/10;
tau_leader = 2*pi/12;
omega_delta = zeros(size(param_space,1),1);
for r = 1:size(param_space,1)
    [rel_phase_r(r),rel_phase(r),rel_phase_error(r),tau(r,:),omega_delta(r)] = ...
        anticipatory_kuramoto(Kvec,epsilon,[tau_1 param_space(r,1)],tau_leader,param_space(r,2),.5,0,0);
end


figure(1)
surf(X,Y,reshape(rel_phase_error,size(X,1),[]),"EdgeColor","black")
xlabel('\tau_2 [rad]')
ylabel('\phi_{leader,0} - \phi_{follower,0} [rad]')
zlabel('\phi')
ylim([0 2*pi])
set(gcf,'color','w')
set(gcf, 'PaperPosition', [0 0 4 3])
set(gcf, 'InvertHardcopy', 'off')
set(gca,'fontsize',14)
if printfigureflag == 1
    f = fullfile(pwd,['rel_phase_final_' char(datetime('now','TimeZone','local','Format','y-MM-d-hhmmss')) '.jpeg']);
    print('-djpeg','-r600',f)
end


figure(2)
surf(X,Y,reshape(omega_delta,size(X,1),[]),"EdgeColor","none")
colormap cool
cb = colorbar;
ylabel(cb, '\Delta \omega [%]')

xlabel('\tau_{follower,2} [degrees]')
ylabel('\phi_{leader,0} - \phi_{follower,0} [degrees]')
xlim([0 2*pi])
set(gca, 'XDir', 'reverse')
set(gca, 'YDir', 'reverse')
set(gca,'XTick',round(0:pi/2:2*pi,2))
set(gca,'XTickLabel',(0:pi/2:2*pi)./pi*180)
set(gca,'YTick',round(0:pi/2:2*pi,2))
set(gca,'YTickLabel',(0:pi/2:2*pi)./pi*180)
ylim([0 2*pi])
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


figure(3)
index = param_space(:,2)==0;
tau_diff = tau_1 + param_space(index,1) - tau_leader;
plot(tau_diff,omega_delta(index,1),'-ok','linewidth',2);
hold on
plot(tau_diff,tau_diff*0,'--k','linewidth',1);
hold off
ylabel('\Delta \omega [%]')
xlabel('\tau_{follower} - \tau_{leader} [degrees]')
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


figure(4)
index = param_space(:,1) == 0;
plot(param_space(index,2),omega_delta(index,1),'-ok','linewidth',2);
hold on
plot(param_space(index,2),param_space(index,2)*0,'--k','linewidth',1);
hold off
ylabel('\Delta \omega [%]')
xlabel('\phi_{leader,0} - \phi_{follower,0} [degrees]')
xlim([0 2*pi])
set(gca,'XTick',round(0:pi/2:2*pi,2))
set(gca,'XTickLabel',(0:pi/2:2*pi)./pi*180)
set(gcf,'color','w')
set(gcf, 'PaperPosition', [0 0 5 4])
set(gcf, 'InvertHardcopy', 'off')
set(gca,'fontsize',14)
if printfigureflag == 1
    f = fullfile(pwd,['omega_drift_by_phi0_' char(datetime('now','TimeZone','local','Format','y-MM-d-hhmmss')) '.jpeg']);
    print('-djpeg','-r600',f)
end
