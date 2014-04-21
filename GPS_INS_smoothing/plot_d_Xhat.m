function r=plot_d_Xhat(Xhat0, Xhat)

K = size(Xhat0);
if K ~= size(Xhat)
    disp('ERROR! Xhat length does not match!');
    return
end
for ii = 1:K
	% log data
	out.t(ii,:) = Xhat0(ii).t;
	out.r_tb_t(ii,:) = (Xhat0(ii).r_tb_t - Xhat(ii).r_tb_t)';
% 	out.v_tb_t(ii,:) = (Xhat0(ii).v_tb_t)';
% 	out.v_tb_b(ii,:) = (Xhat0(ii).R_b2t'*Xhat0(ii).v_tb_t)';
%     out.R_b2t(:,:,ii)  = Xhat0(ii).R_b2t;
% 	out.ba_b(ii,:) = Xhat0(ii).ba_b';
% 	out.bg_b(ii,:) = Xhat0(ii).bg_b';
% 	%out.Pxx(ii,:) = diag(nav.dx.Pxx_minus);
% 	xt = Xhat0(ii).R_b2t*[1;0;0];
% 	out.alphahat(ii,:) = atan2(xt(2),xt(1));
end

r = std(out.r_tb_t);

h = figure(2);
color = get(gca,'ColorOrder');
subplot(3,1,1);
plot(out.t, out.r_tb_t(:,1),'--*b');
title(['NED position'] );
hold on;
grid on;
ylabel('x (m)');
subplot(3,1,2);
plot(out.t, out.r_tb_t(:,2),'--*b');
hold on;
grid on;
ylabel('y (m)');
subplot(3,1,3);
plot(out.t, out.r_tb_t(:,3),'--*b');
hold on;
grid on;
ylabel('z (m)');