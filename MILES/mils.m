function [Xhat,Zhat] = mils(A,B,y,p)
%
% [Xhat,Zhat] = mils(A,B,y,p) produces p pairs of optimal solutions to 
%               the mixed integer least squares problem min_{x,z}||y-Ax-Bz||, 
%               where x and z are real and integer vectors, respectively.
%
% Input arguments:
%    A ---- m by k real matrix
%    B ---- m by n real matrix
%          [A,B] has full column rank
%    y ---- m-dimensional real vector
%    p ---- the number of optimal solutions and its default value is 1
%
% Output arguments:
%    Xhat - k by p real matrix
%    Zhat - n by p integer matrix (in double precision). 
%           The pair {Xhat(:,j),Zhat(:,j)} is the j-th optimal solution
%           i.e., its residual is the j-th smallest, so
%           ||y-A*Xhat(:,1)-B*Zhat(:,1)||<=...<=||y-A*Xhat(:,p)-B*Zhat(:,p)||

% Copyright (c) 2006. Xiao-Wen Chang and Tianyang Zhou
% Version 1.0, October 2006.


% Check input arguments
if nargin < 3 % input error
    error('Not enough input arguments!')
end

if nargin <= 3
    p = 1;
end

if p <= 0 % input error
    error('Third input argument must be an integer bigger than 0!')
end

[m,k] = size(A);
[m2,n] = size(B);
if m ~= m2 || m ~= size(y,1) || size(y,2) ~= 1  % input error
    error('Input arguments have a matrix dimension error!')
end

if rank([A,B]) < k+n
    error('The pair of the two matrices is rank defficient!')
end

[Q,R] = qr(A);
Q_A = Q(:,1:k); 
Q_Abar = Q(:,k+1:m);
R_A = R(1:k,:);

% Compute the p optimal integer least squares solutions
Zhat = ils(Q_Abar'*B,Q_Abar'*y, p);

% Compute the corresponding real least squares solutions
Xhat  = R_A\(Q_A'*(y*ones(1,p)-B*Zhat));
