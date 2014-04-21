function plot_Xhat(nav, Xhat_k, K, c_num)
r2d = 180/pi;
for ii = 1:K
	% log data
	out.t(ii,:) = Xhat_k(ii).t;
	out.r_tb_t(ii,:) = (Xhat_k(ii).r_tb_t)';
	out.v_tb_t(ii,:) = (Xhat_k(ii).v_tb_t)';
	out.v_tb_b(ii,:) = (Xhat_k(ii).R_b2t'*Xhat_k(ii).v_tb_t)';
    out.R_b2t(:,:,ii)  = Xhat_k(ii).R_b2t;
	out.ba_b(ii,:) = Xhat_k(ii).ba_b';
	out.bg_b(ii,:) = Xhat_k(ii).bg_b';
	%out.Pxx(ii,:) = diag(nav.dx.Pxx_minus);
	xt = Xhat_k(ii).R_b2t*[1;0;0];
	out.alphahat(ii,:) = atan2(xt(2),xt(1));
end

h = figure(1);
color = get(gca,'ColorOrder');
subplot(3,1,1);
plot(out.t, out.r_tb_t(:,1)-nav.param.base_pos(1),'--*b', 'Color', color(c_num,:));
title(['NED position'] );
hold on;
grid on;
ylabel('x (m)');
subplot(3,1,2);
plot(out.t, out.r_tb_t(:,2)-nav.param.base_pos(2),'--*b', 'Color', color(c_num,:));
hold on;
grid on;
ylabel('y (m)');
subplot(3,1,3);
plot(out.t, out.r_tb_t(:,3)-nav.param.base_pos(3),'--*b', 'Color', color(c_num,:));
hold on;
grid on;
ylabel('z (m)');

h = figure(2);
color = get(gca,'ColorOrder');
subplot(3,1,1);
plot(out.t, out.v_tb_t(:,1),'--*b', 'Color', color(c_num,:));
title(['NED velocity'] );
hold on;
grid on;
ylabel('vx (m)');
subplot(3,1,2);
plot(out.t, out.v_tb_t(:,2),'--*b', 'Color', color(c_num,:));
hold on;
grid on;
ylabel('vy (m)');
subplot(3,1,3);
plot(out.t, out.v_tb_t(:,3),'--*b', 'Color', color(c_num,:));
hold on;
grid on;
ylabel('vz (m)');

h = figure(3);
color = get(gca,'ColorOrder');
plot(out.t, limit_pi(out.alphahat)*r2d,'--*b', 'Color', color(c_num,:));
title(['Yaw'] );
hold on;
grid on;

h = figure(4);
color = get(gca,'ColorOrder');
plot(out.r_tb_t(:,2)-nav.param.base_pos(2), out.r_tb_t(:,1)-nav.param.base_pos(1), '--*b', 'Color', color(c_num,:));
title(['Tangent Plane Trajectory'] );
ylabel('North (m)');
xlabel('East (m)');
hold on;
grid on;