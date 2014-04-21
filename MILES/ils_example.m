%
% A small example to run function ils.m
%

% Construct data
m = 5;
n = 3; 
B = randn(m,n);
z_true = [1; -2; 3];
y = B*z_true + 1.e-3*randn(m,1);
p = 2;

% Find p optimal solution to the ILS problem min_{z}||y-Bz|| 
Zhat = ils(B,y,p);

display('The true integer parameter vector')
z_true

display('The two integer least squares estimates')
Zhat

