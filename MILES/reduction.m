function [R,Z,y] = reduction(B,y)
%
% [R,Z,y] = reduction(B,y) computes the LLL-QRZ factorization:
%           Q'*B*Z = [R; 0] and computes Q'*y. The orthogonal matrix Q 
%           is not produced. Its goal is to reduce a general integer
%           least squares problem to an upper triangular one.
%
% Input arguments:
%    B --- m by n real matrix with full column rank
%    y --- m-dimensional real vector to be transformed to Q'y
%
% Output arguments:
%    R --- n by n LLL-reduced upper triangular matrix
%    Z --- n by n unimodular matrix, i.e., an integer matrix with |det(Z)|=1
%    y --- m-vector transformed from the input y by Q', i.e., y := Q'*y

% Copyright (c) 2006-2012. Xiao-Wen Chang, Xiaohu Xie and Tianyang Zhou
% Version 2.0, October 2011.



% Check input arguments
if nargin < 2  % input error
    error('Not enough arguments!')
end

[m,n] = size(B);
if m < n  % input error
    error('Input matrix is column-rank deficient!')
end

[m2,n2] = size(y);
if m ~= m2 | n2 ~= 1  % input error
    error('Input arguments have a dimension error!')
end

d2=1;

% QR with minimum-column pivoting
[R,piv,y] = qrmcp(B,y);

% Obtain the permutation matrix Z
Z = zeros(n,n);
for j = 1:n
    Z(piv(j),j) = 1;
end


% Perform partial LLL reduction on R

k = 2;

while k <= n
    
    k1 = k-1;
    zeta = round(R(k1,k)/R(k1,k1));  
    alpha = R(k1,k)-zeta*R(k1,k1);  

    if R(k1,k1)^2 > (1+1.e-10)*(alpha^2+R(k,k)^2)   
        if zeta ~= 0
            % Perform a size reduction on R(k-1,k)
            R(k1,k) = alpha;
            R(1:k-2,k) = R(1:k-2,k) - zeta*R(1:k-2,k-1);
            Z(:,k) = Z(:,k) - zeta*Z(:,k-1);  
            
            % Perform size reductions on R(1:k-2,k)
            for i = k-2:-1:1
                zeta = round(R(i,k)/R(i,i));  
                if zeta ~= 0
                    R(1:i,k) = R(1:i,k) - zeta*R(1:i,i);  
                    Z(:,k) = Z(:,k) - zeta*Z(:,i);  
                end
            end
        end
        
        % Permute columns k-1 and k of R
        R(1:k,[k1,k]) = R(1:k,[k,k1]);
        Z(:,[k1,k]) = Z(:,[k,k1]);
        
        % Bring R back to an upper triangular matrix by a Givens rotation
        r = sqrt(R(k1,k1)^2+R(k,k1)^2);
        c = R(k1,k1)/r;
        s = R(k,k1)/r;
        G = [c,s; -s,c];
        R(k1,k1) = r;
        R(k,k1) = 0;
        R([k1,k],k:n) = G*R([k1,k],k:n);
        
        % Apply the Givens rotation to y
        y([k1,k]) = G*y([k1,k]);
        
        if k > 2
            k = k-1;
        end
        
    else
        
        k=k+1;
        
    end
end



