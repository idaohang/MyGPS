function [R,piv,y] = qrmcp(B,y)
%
% [R,piv,y] = qrmcp(B,y) computes the QR factorization of B with 
%             minimum-column pivoting: Q'BP = [R; 0] and computes Q'*y. 
%             The orthogonal matrix Q is not produced. 
%
% Input arguments:
%    B --- m by n real matrix to be factorized
%    y --- m-dimensional real vector to be transformed to Q'y
%
% Output arguments:
%    R --- n by n real upper triangular matrix
%    piv - n-vector storing the information of the permutation matrix P
%    y --- m-vector transformed from the input y by Q, i.e., y := Q'*y 

% Copyright (c) 2006-2011. Xiao-Wen Chang, Xiaohu Xie and Tianyang Zhou
% Version 2.1, July 2012


% Check input arguments
[m,n] = size(B);
if m < n  % input error
    error('Matrix to be factorized is column-rank deficient!')
end;

[m2,n2] = size(y);
if m ~= m2 | n2 ~= 1  % input error
    error('Input arguments have a dimension error!')
end

% Initialization
colnormB = zeros(2,n);
piv = 1:n;

% Compute the 2-norm squared of each column of B 
for j = 1:n
    colnormB(1,j) = (norm(B(:,j)))^2;
end


for k = 1:n

    % Find the column with minimum 2-norm in B(k+1:m,k+1:n)
    [minnorm, i] = min(colnormB(1,k:n)-colnormB(2,k:n));
    q = i + k - 1;
    
    % Column interchange
    if q > k
        piv([k,q]) = piv([q,k]);
        colnormB(:,[k,q]) = colnormB(:,[q,k]);
        B(:,[k,q]) = B(:,[q,k]);
    end

    % Compute and apply the Householder transformation  I-tau*v*v'
    if norm(B(k+1:m,k)) > 0, % otherwise no Householder transformation is needed
	    v = B(k:m,k);
        rho = norm(v);     
	    if (v(1) >= 0) 
            rho = -rho;
        end
        v(1) = v(1) - rho;   % B(k,k)+sgn(B(k,k))*norm(B(k:n,k))
        tao = -1/(rho*v(1));
        B(k,k) = rho;
        B(k:m,k+1:n) = B(k:m,k+1:n) - tao*v*(v'*B(k:m,k+1:n));
    
        % Update y by the Householder transformation
        y(k:m) = y(k:m,:) - tao*v*(v'*y(k:m));
    end
  
    % Update colnormB(2,k+1:n)
    colnormB(2,k+1:n) = colnormB(2,k+1:n) + B(k,k+1:n).*B(k,k+1:n);

end

R = triu(B(1:n,1:n));


