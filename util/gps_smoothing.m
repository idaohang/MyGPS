function [out, Xhat_k] = gps_smoothing(nav,data)
% dimension of states
N = nav.param.X_STATES;

% number of states
K = nav.buff.kdim+1;

% dimension of correction
l = N*K;

% initial conditions
Xhat_c = [nav.buff.xhat_k0;nav.buff.xhat_k];
Xhat_t = Xhat_c; % initialize temp Xhat_t

itc = 1;
maxit =20;
tol = 0.1;
alpha = 1.d-4;
[rc, Jc, gc] = gps_smth_residual_jacobian(nav, Xhat_c, data);
n_rc0 = norm(rc)
while(norm(gc)>tol && itc<=maxit)
    lambda = 1;
    % Gauss-Newton
    [Q, R] = qr(Jc);
    de = -Q'*rc;
    d = de(1:l, :);
    R = R(1:l, :);
    dc = R\d;
    dgn = lambda*dc;
    %n_dgn = norm(dgn);
    delta = reshape(dgn, N, K)';
    for k = 1:size(Xhat_c,1)
        Xhat_t(k) = gpsins_correct_xhat(nav.param, Xhat_c(k), delta(k,:)');
    end
    
    [rt, Jt, gt] = gps_smth_residual_jacobian(nav, Xhat_t, data);
    iarm=0;
    itc=itc+1;
    %
    % Goal for sufficient decrease
    %    
    rgoal = norm(rc) - alpha*lambda*(gc'*dc); % dgn=lambda*dc
    Armijo = 0;
    while(norm(rt)>rgoal)        
        iarm=iarm+1;
        lambda=lambda/5;
        rgoal = norm(rc) - alpha*lambda*(gc'*dc);
        dgn = lambda*dc;
        delta = reshape(dgn, N, K)';
        for k = 1:size(Xhat_c,1)
            Xhat_t(k) = gpsins_correct_xhat(nav.param, Xhat_c(k), delta(k,:)');
        end
        [rt, Jt, gt] = gps_smth_residual_jacobian(nav, Xhat_t, data);
        if(iarm > 10)
            disp(' Armijo error in Gauss-Newton')
            Armijo = 1;
            break;
        end
    end
    if Armijo
        break;
    else
        Xhat_c = Xhat_t;
        [rc, Jc, gc] = gps_smth_residual_jacobian(nav, Xhat_c, data);
        norm(rc)
    end
end
Xhat_k = Xhat_c;

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
subplot(3,1,1);
plot(out.t, out.r_tb_t(:,1)-nav.param.true_pos(1),'*r');
title(['NED position'] );
hold on;
grid on;
ylabel('x (m)');
subplot(3,1,2);
plot(out.t, out.r_tb_t(:,2)-nav.param.true_pos(2),'r');
hold on;
grid on;
ylabel('y (m)');
subplot(3,1,3);
plot(out.t, out.r_tb_t(:,3)-nav.param.true_pos(3),'r');
hold on;
grid on;
ylabel('z (m)');

end