function [error, variance, stdev] = check_result_error(out, state_out_k, shift)
if nargin <3
    shift = 0;
end

while state_out_k(1,1) < out(1,1)
    state_out_k(1,:) = [];
end

[m,n] = size(state_out_k);
error = zeros(m,3);

if shift && norm(out(1,2:3) - state_out_k(1,2:3))>0.05
    dis = zeros(length(out), 1);
    
    for j = 1:length(out)
        dis(j,:) = norm(out(j,2:3) - state_out_k(1,2:3));
    end
    
    [ind, ind] = min(dis);
    tm_shift = state_out_k(1,1) - out(ind,1)
    if abs(tm_shift)<=0.01
        out(:,1) = out(:,1) + tm_shift;
    end
end

for i = 1:m
    tm = state_out_k(i,1);
    [ind, ind] = min(abs(out(:,1)-tm));
    error(i,:) = out(ind, 2:4) - state_out_k(i,2:4);
end
variance = var(error);
stdev = std(error);
error = mean(error);

