function Zhat = ils(B,y,p)
%
% Zhat = ils(B,y,p) produces p optimal solutions to the integer least 
%        squares problem min_{z}||y-Bz||
%
% Input arguments:
%    B ---- m by n real matrix with full column rank
%    y ---- m-dimensional real vector
%    p ---- the number of optimal solutions and its default value is 1
%
% Output arguments:
%    Zhat - n by p integer matrix (in double precision). Its j-th column 
%           is the j-th optimal solution, i.e., its residual is the j-th
%           smallest, so ||y-B*Zhat(:,1)|| <= ...<= ||y-B*Zhat(:,p)||

% Copyright (c) 2006. Xiao-Wen Chang and Tianyang Zhou 
% Version 1.0, October 2006.
 

% Check input arguments
if nargin < 2 % input error
    error('Not enough input arguments!')
end

if nargin < 3
    p = 1;
end

if p <= 0 % input error
    error('Third input argument must be an integer bigger than 0!')
end

[m,n] = size(B);
if m < n || m ~= size(y,1) || size(y,2) ~= 1  % input error
    error('Input arguments have a matrix dimension error!')
end

if rank(B) < n
	error('Matrix is rank defficient!')
end

% Reduction
[R,Z,y] = reduction(B,y);

% Search
Zhat = search(R,y(1:n),p);

% Perform the unimodual transformation to obtain the optimal solutions
Zhat = Z*Zhat;
