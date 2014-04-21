% y = Ax + Bz + n

epoch_num = 3;
state_num = 3;
sat_num   = 6;
% A = zeros( epoch_num*sat_num, state_num*epoch_num);
B = zeros( epoch_num*sat_num, state_num)
A1 = [1,1,1,0,0,0,0,0,0; 0,0,0,1,1,1,0,0,0;0,0,0,0,0,0,1,1,1 ];
A2 = A1*2;
A3 = A1*3;
A4 = A1*4;
A5 = A1*5;
A6 = A1*6;

A = [A1; A2; A3; A4; A5; A6];
for i = 1:sat_num
    B(i*3-2:i*3, i) = ones(3,1);
end

[Q, ~] = qr(B);
QQ = Q(:,sat_num+1:end);

AA = QQ'*A;


[Qs, ~] = qr(ones(3,1));
Qs = Qs(:,2:end);

As = [Qs'*A1; Qs'*A2; Qs'*A3; Qs'*A4; Qs'*A5; Qs'*A6];