function nav = msc_gps_buff_init(nav, K)
nav.Xhat(K,1) = nav.xhat;
nav.Dx = zeros(K*nav.param.X_STATES,1);
nav.Pxx = [];
end