% We are given a covariance matrix, and we’d like to generate n vectors, whose coordinates’ covariance matrix is equal to given. We also like each coordinate to have given mean.
% n – number of vectors to generate
% covMatrix – given covariance matrix
% offset – vector of means of each coordinate
function numbers = randomcov(n, covMatrix, offset)
    x=randn(n,length(covMatrix(1,:)));
    numbers = x*inv(chol(cov(x)))*chol(covMatrix);
    for i=1:length(offset)
        numbers(:,i) = offset(i)+numbers(:,i)-mean(numbers(:,i));
    end
end