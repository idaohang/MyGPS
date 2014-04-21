
function [r, J, g, JN, Ndim] = gps_smth_res_J_shift(nav, data, Xhat_k, Phase_meas)
global Dual_Freq;
dual_freq = Dual_Freq;
p = nav.param;  % p.X_POS is the index for tagent position in std state vector

% dimensions
K = nav.buff.kdim;     % number of xhat
Ndim = Phase_meas.num_ph_sv - 1;

% initialize residual r
SigmaP = eye(3); %chol(inv(nav.buff.Pxx_k0(p.X_POS, p.X_POS)));
temp = gpsins_compute_dx(nav.param, Xhat_k(1), Xhat_k(1)); % the last input changed from nav.xhat0
rP = SigmaP*temp(p.X_POS);

% initialize Jacobian J(Xhat)
JP = [ SigmaP ];
JYs = zeros( Ndim, 3 );
JN  = zeros( Ndim, 1 ); % for the Jacobian of N
rYs = zeros(Ndim, 1);
Rs = zeros(Ndim, 1);
JYs_wd = zeros( Ndim, 3 );
JN_wd  = zeros( Ndim, 1 ); % for the Jacobian of N
rYs_wd = zeros(Ndim, 1);
Rs_wd = zeros(Ndim, 1);

Pos_ecef = zeros(K*3,1);
for n=1:K
    Pos_ecef(3*n-2:3*n) = nav.param.ecef_p_b + nav.param.R_ned2ecef*Xhat_k(n+1).r_tb_t;
end

for i = 1:Ndim
    [H_abs, R_abs, dy_abs, J_n] = abs_phase_output(nav.param, Pos_ecef, data.gps(1:end), i+1, K, Phase_meas);
    JYs(i,:) = H_abs(1,:);
    JN(i) = J_n(1,:);
    rYs(i) = dy_abs(1,:);
    Rs(i) = R_abs(1,:);
    if dual_freq
        JYs_wd(i,:) = H_abs(2,:);
        JN_wd(i) = J_n(2,:);
        rYs_wd(i) = dy_abs(2,:);
        Rs_wd(i) = R_abs(2,:);
    end
end
if dual_freq
    JYs = [JYs; JYs_wd];
    JN = [JN; JN_wd];
    rYs = [rYs; rYs_wd];
    Rs = [Rs; Rs_wd];
end

SigmaR = chol( diag(1./Rs) );
JYs = SigmaR*JYs;
rYs = SigmaR*rYs;
JN  = [zeros(3,Ndim*(1+dual_freq));diag(SigmaR*JN)];

% output
r = [rP; rYs];
J = [JP; JYs];

g = J'*r;

%% No marginalization
% Ns = size(p.X_POS,2);    % size of trajectory shift
% Ys = Phase_meas.num_ph_meas;
%
% rYs = zeros(Ys, 1);
%
% % initialize Jacobian J(Xhat)
% JP = [SigmaP zeros( Ns, Ndim )];
% JYs = zeros( Ys, Ns+Ndim );
%
% y_s = 1;
% for k = 1:size(nav.buff.index,1)
%     xhat_k = Xhat_k(k+1);
%
%     ii_gps = nav.buff.index(k,1);
%     % sample GPS measurements
%     pos_ecef = nav.param.ecef_p_b + nav.param.R_ned2ecef*xhat_k.r_tb_t;
%     [Hs, Rs, dys] = gpsins_gps_output_int(nav.param, pos_ecef, data.gps(ii_gps), Phase_meas);
%
%     [rYs, JYs, y_s] = accumulateRNs(rYs, JYs, Hs, Rs, dys, y_s, k, Ndim, p.wave_l1); %k
% end
%
% % check dimensions
% if k ~= K
%     error('invalid k-steps')
% end
%
% % output
% r = [rP; rYs];
% J = [JP; JYs];
% JN = J(:,4:end);
% J  = J(:,1:3);
% g = J'*r;

function [rY, JY, y] = accumulateRs(rY, JY, H, R, dy, y, k)
N = size(H, 2);
m = size(dy, 1);
if ~isempty(dy)
    SigmaR = chol(inv(R));
    rY(y:(y+m-1), 1) = -SigmaR*dy;
    JY(y:(y+m-1), 1:3) = SigmaR*H(:,7:9);
    y = y + m;
end

function [rY, JY, y] = accumulateRNs(rY, JY, H, R, dy, y, k, Ndim, lambda)
N = size(H, 2);
m = size(dy, 1);
if ~isempty(dy)
    SigmaR = chol(inv(R));
    rY(y:(y+m-1), 1) = SigmaR*dy;
    JY(y:(y+m-1), 1:3) = SigmaR*H(:,7:9);
    JY(y:(y+m-1), (end-Ndim+1):end) = SigmaR*lambda*eye(Ndim);
    y = y + m;
end
