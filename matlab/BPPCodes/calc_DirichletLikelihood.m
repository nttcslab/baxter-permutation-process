function result = calc_DirichletLikelihood(EachClusterCount, alpha)

alpha = alpha * ones(size(EachClusterCount,1),1);
result = 0;

for ii=1:size(EachClusterCount,2)
    if sum(EachClusterCount(:,ii),1)~=0
        result = result + calc_marginalDir(EachClusterCount(:,ii), alpha);
    end
end





function log_marginalDir = calc_marginalDir(N, alpha)


log_marginalDir = gammaln(sum(alpha,1)) + sum(gammaln(N+alpha),1) - gammaln(sum(N,1)+sum(alpha,1)) - sum(gammaln(alpha),1);






