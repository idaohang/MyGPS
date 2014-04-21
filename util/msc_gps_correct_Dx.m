function [nav] = msc_gps_correct_Dx(nav, H, R, dy, buff_length)
if ~isempty(dy)
	%I = eye(size(dx.Pxx_minus));
    N = nav.param.X_STATES;
    I = eye(N*buff_length);
    
    nav.Pxx = [];
    for n=1:buff_length
        nav.Pxx = blkdiag(nav.Pxx,nav.dx_buff(n).Pxx_minus);
    end
    
	HPHtR = H*nav.Pxx*H' + R;
	K = nav.Pxx*H'/(HPHtR);
	nav.Dx = K*(dy - H*nav.Dx);
	nav.Pxx = (I - K*H)*nav.Pxx*(I - K*H)' + K*R*K';
    
%     for n = 1:buff_length
%         nav.dx_buff(n).Pxx_minus = nav.Pxx(N*(n-1)+1:N*n,N*(n-1)+1:N*n);
%     end   
    nav.dx.Pxx_minus = nav.Pxx(N*(buff_length-1)+1:N*buff_length,N*(buff_length-1)+1:N*buff_length);    
end
