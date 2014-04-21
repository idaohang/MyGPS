function []=gps_plot_Xhat(p, Xhat, dx_buff, K, fig_number)
for ii = 1:K
    % log data
    tout.t(ii,:) = Xhat(ii).t;
    tout.r_tb_t(ii,:) = (Xhat(ii).r_tb_t)';
    tout.v_tb_t(ii,:) = (Xhat(ii).v_tb_t)';
    tout.v_tb_b(ii,:) = (Xhat(ii).R_b2t'*Xhat(ii).v_tb_t)';
    tout.R_b2t(:,:,ii)  = Xhat(ii).R_b2t;
    tout.ba_b(ii,:) = Xhat(ii).ba_b';
    tout.bg_b(ii,:) = Xhat(ii).bg_b';
    tout.Pxx(ii,:) = diag(dx_buff(ii).Pxx_minus);
    xt = Xhat(ii).R_b2t*[1;0;0];
    tout.alphahat(ii,:) = atan2(xt(2),xt(1));
end

h = figure(fig_number);
subplot(3,1,1);
plot(tout.t, tout.r_tb_t(:,1)-p.true_pos(1), '*', tout.t,[3*sqrt(tout.Pxx(:,7)),-3*sqrt(tout.Pxx(:,7))], '--b');
title(['NED position'] );
grid on;
ylabel('x (m)');
subplot(3,1,2);
plot(tout.t, tout.r_tb_t(:,2)-p.true_pos(2), '*', tout.t,[3*sqrt(tout.Pxx(:,8)),-3*sqrt(tout.Pxx(:,8))], '--b');
grid on;
ylabel('y (m)');
subplot(3,1,3);
plot(tout.t, tout.r_tb_t(:,2)-p.true_pos(2), '*', tout.t,[3*sqrt(tout.Pxx(:,9)),-3*sqrt(tout.Pxx(:,9))], '--b');
grid on;
ylabel('z (m)');