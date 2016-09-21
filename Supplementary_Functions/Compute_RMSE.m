function RMSE = Compute_RMSE(M,M_est)

[m, n] = size(M);
RMSE = sqrt(sum(sum((M - M_est).^2))/(m*n));
