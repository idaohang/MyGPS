% This script is written by Yiming Chen to show Sheng's idea of first
% marginalizing B matrix will not work for MILS, since 'z is an integer
% vector' is a stronger constraint.s

clear all;
clc;

%
% A small example to run function mils.m
%

% Construct data
m = 7;
k = 2;
n = 3;
p = 3;

A = randn(m,k)
B = randn(m,n)
y = randn(m,1)

% To minimize \| y - Ax - Bz \|

display('Three pairs of optimal least squares solutions')
[X,Z] = mils(A,B,y,p)

rk = rank(B)
[Qb,Rb]=qr(B)
Q2 = Qb(:,rk+1:end)
% check Q2'*B == 0
if norm(Q2'*B) > 1e-14
    error('Q2 does not marginalize B matrix')
end
new_y = Q2' * y
new_H = Q2' * A
new_x = new_H \ new_y

display('The estimation error is')
error = norm( new_x - X(:,1) )