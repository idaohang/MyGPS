
function [r, J, g, gps_res] = gps_smth_res_J_mils(nav, data, Xhat_k, Phase_meas)
global Dual_Freq;
dual_freq = Dual_Freq;
p = nav.param;

% dimensions
Ns = nav.param.X_STATES;
K = nav.buff.kdim;     % number of xhat
Y = nav.buff.ydim;     % number of measurements, L1 code
YN = Phase_meas.num_ph_meas;
Ndim = (Phase_meas.num_ph_sv - 1);

% initialize residual r(Xhat)
SigmaP = chol(inv(nav.buff.Pxx_k0));
rP = SigmaP*gpsins_compute_dx(nav.param, Xhat_k(1), nav.xhat0);


rQ = zeros(Ns*K, 1);
rY = zeros(Y, 1);
rYN = zeros(YN, 1);

%fprintf(1, 'P_0(7,7)=%-12g P_0(8,8)=%-12g P_0(9,9)=%g\n', Pxx_k0(7,7), Pxx_k0(8,8), Pxx_k0(9,9));

% initialize Jacobian J(Xhat)
JP = [SigmaP zeros( Ns, K*Ns )];
JQ = zeros( Ns*K, Ns*(K+1));
JY = zeros( Y, Ns*(K+1) );
JYN = zeros( YN, Ns*(K+1)+Ndim*(dual_freq+1) );
%JYN = zeros( Y, Ns*(K+1) );

gps_res.t = zeros(size(nav.buff.index,1),1);
gps_res.r = zeros(size(nav.buff.index,1),Ndim);
gps_res.prn_list = Phase_meas.prn_phase;
y = 1;
y_N = 1;
for k = 1:size(nav.buff.index,1) 
    xhat_k_1 = Xhat_k(k);
    xhat_k = Xhat_k(k+1);
    
    ii_gps = nav.buff.index(k,1);
    % sample GPS measurements
    pos_ecef = nav.param.ecef_p_b + nav.param.R_ned2ecef*xhat_k.r_tb_t;
    [H, R, dy] = gpsins_gps_output(nav.param, pos_ecef, data.gps(ii_gps), 1, dual_freq);
    [HN, RN, dyN] = gpsins_gps_output_int(nav.param, pos_ecef, data.gps(ii_gps), Phase_meas);
    
    [rY, JY, y] = accumulateR(rY, JY, H, R, dy, y, k); %k   
    [rYN, JYN, y_N] = accumulateRN(rYN, JYN, HN, RN, dyN, y_N, k, Ndim, p); %k

    gps_res.t(k,1) = xhat_k.t;
    gps_res.r(k,:) = (-dyN(1:Ndim))';
    
    % repropagate and accumulate Q
    [xbar_k, Phi_k, Qd_k, nav.buff.Qd_k0] = gps_smth_iterate(nav.param, xhat_k_1, nav.buff.Qd_i0, [nav.buff.index(k,2),nav.buff.index(k,3)], data);
    [rQ, JQ] = accumulateQ(nav.param, rQ, JQ, xbar_k, xhat_k, Phi_k, Qd_k, k);
    
    %fprintf(1, '\nQ%g(7,7)=%-12g, Q%g(8,8)=%-12g, Q%g(9,9)=%g\n', k, Qd_k(7,7), k,Qd_k(8,8), k,Qd_k(9,9));
end

% check dimensions
if k ~= K
    error('invalid k-steps')
end
% if y-1 ~= Y
%     error('invalid measurement dimensions');
% end

% output
JP = [JP, zeros( Ns, Ndim*(dual_freq+1) )];
JQ = [JQ, zeros( Ns*K, Ndim*(dual_freq+1))];
JY = [JY, zeros( Y, Ndim*(dual_freq+1))];
r = [rP*0; rQ; rY; rYN];
J = [JP; JQ; JY; JYN];
g = J'*r;
gps_res = [gps_res.t, gps_res.r];


function [rQ, JQ] = accumulateQ(p, rQ, JQ, xbar_k, xhat_k, Phi_k, Qd_k, k)
N = size(Phi_k, 1);
SigmaQ = chol(inv(Qd_k));
a_k = -gpsins_compute_dx(p, xhat_k, xbar_k);
rQ(((k-1)*N+1):(k*N), 1) = SigmaQ*a_k;
JQ(((k-1)*N+1):(k*N), ((k-1)*N+1):(k*N)) = SigmaQ*Phi_k;
JQ(((k-1)*N+1):(k*N), (k*N+1):((k+1)*N)) = -SigmaQ;


function [rY, JY, y] = accumulateR(rY, JY, H, R, dy, y, k)
N = size(H, 2);
m = size(dy, 1);
if ~isempty(dy)
    SigmaR = chol(inv(R));
    rY(y:(y+m-1), 1) = -SigmaR*dy;
    JY(y:(y+m-1), (k*N+1):((k+1)*N)) = SigmaR*H;
    y = y + m;
end

function [rY, JY, y] = accumulateRN(rY, JY, H, R, dy, y, k, Ndim, p)
global Dual_Freq;
dual_freq = Dual_Freq;
N = size(H, 2);
m = size(dy, 1);
if ~isempty(dy)
    SigmaR = chol(inv(R));
    rY(y:(y+m-1), 1) = SigmaR*dy;
    JY(y:(y+m-1), (k*N+1):((k+1)*N)) = -SigmaR*H;
    JY(y:(y+m/(dual_freq+1)-1), (end-Ndim*(dual_freq+1)+1):end-Ndim*dual_freq) = SigmaR(1:Ndim,1:Ndim)*p.wave_l1*eye(Ndim);
    if dual_freq
        JY(y+m/(dual_freq+1):(y+m-1), (end-Ndim+1):end) = SigmaR(Ndim+1:end,Ndim+1:end)*p.wave_wd*eye(Ndim);
    end
    y = y + m;
end
